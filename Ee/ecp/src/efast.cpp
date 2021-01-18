#include<Rcpp.h>
#include<math.h>

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

// user includes
#include<algorithm>
#include<iostream>
#include<vector>
#include<set>
// #include<random>

using namespace std;
using namespace Rcpp;

// declarations
extern "C"{
SEXP eFastC_delta(SEXP Z_, SEXP K_, SEXP delta_, SEXP alpha_, SEXP verbose_);
SEXP eFastC(SEXP Z_, SEXP K_, SEXP minsize_, SEXP alpha_, SEXP verbose_);
SEXP ksFastC_delta(SEXP Z_, SEXP K_, SEXP minsize_, SEXP verbose_);
SEXP ksFastC(SEXP Z_, SEXP K_, SEXP minsize_, SEXP verbose_);
}

double dst(const NumericVector& X, const NumericVector& Y, double alpha);
std::vector<double> MEAN_VAR(const std::vector<double>& X);
double delta_sum(NumericMatrix& X, int a, int b, double alpha);
double dist_X(NumericMatrix& X, double alpha);
double dist_XY(NumericMatrix& X, NumericMatrix&Y, double alpha);
std::vector<std::vector<int> > find_locations(const NumericMatrix& A);
double dist_ks(std::vector<double>& X, std::vector<double>& Y);

// definitions

SEXP eFastC_delta( SEXP Z_, SEXP K_, SEXP delta_, SEXP alpha_, SEXP verbose_){
BEGIN_RCPP
	
	//Convert SEXP variables to types Rcpp/C++ knows
	int K = as<int>(K_), delta = as<int>(delta_);
	double alpha = as<double>(alpha_);
	bool verbose = as<bool>(verbose_);
	NumericMatrix Z(Z_);
	
	int N = Z.nrow(); //Get number of observations

	// Reset K to be maximum number of possible segments
	if(N/(delta+1) < K)
		K = N/(delta+1);

	//Since the update of the goodness-of-fit value with k change points only depends
	//upon the goodness-of-fit values with (k-1) change points we can use a 2xN matrix.
	//However, we will need a KxN matrix to store the change point locaitons.
	NumericMatrix FF(2,N), A(K,N); //Create matrices used to store goodness-of-fit values and 
	//change point locations

	//Goodness-of-fit values for when outer loop is at the end of the time series.
	//GOF[k] corresponds to the goodness-of-fit value  for the entire time series when
	//using k+1 change points.
	NumericVector GOF(K); 

	std::fill(FF.begin(), FF.end(), R_NegInf); //Fill FF matrix with negative infinity
	std::fill(A.begin(), A.end(), -1); //Fill A matrix with -1 because of C arrays being 0 indexed

	//Counters and distances

	//Sum of within sample for left and right segments as well as between segment distances
	double dll = 0.0, drr = 0.0, dlr = 0.0;
	//Counter for number of terms in each sum
	int cll = 0, crr = 0, clr = 0;

	//Precalculations
	if(verbose)
		Rcout<<"Starting pre-processing"<<std::endl;
	int minsize = delta + 1;
	
	//Calculations for complete portion of statistics in delta neighborhoods
	NumericVector DLL(N), DRR(N), DLR(N);
	//DRR[s] = within sum for Z[s+1], Z[s+2], ..., Z[s+delta]
	//DLL[s] = within sum for Z[s-delta+1], Z[s-delta+2], ..., Z[s]
	//DLR[s] = between sum for the sets used to calculate DLL[s] and DRR[s]
	for(int s = delta; s < N; ++s){
		DLL[s] = delta_sum(Z, s-delta+1, s, alpha);
		if(s >= delta)//avoid array out of bounds error
			DRR[s-delta] = DLL[s];
	}

	//Calculate DLR array in O(delta*N) time

	NumericMatrix Left(N,2), Right(N,2);
	//Left(i,0) = sum of distances of Z[i] to {Z[i-1], Z[i-2], ..., Z[i-delta]}
	//Right(i,0) = sum of distances of Z[i] to {Z[i+1], Z[i+2], ..., Z[i+delta]}
	//Left(i,1) = sum of distances of Z[i] to {Z[i-delta], Z[i-delta-1], ..., Z[i-2*delta+1]}
	//Right(i,1) = sum of distances of Z[i] to {Z[i+delta], Z[i+delta+1], ..., Z[i+2*delta-1]}

	for(int i = delta; i < N-delta; ++i){
		for(int j1 = i-delta; j1 < i; ++j1)
			Left(i,0) += dst(Z(i,_), Z(j1,_), alpha);
		for(int j2 = i+1; j2 < i+delta+1; ++j2)
			Right(i,0) += dst(Z(i,_), Z(j2,_), alpha);

		//Avoid array out of bounds errors
		if(i >= 2*delta-1)
			for(int j3=i-2*delta+1; j3<i-delta+1; ++j3)
				Left(i,1) += dst(Z(i,_), Z(j3,_), alpha);
		if(i+2*delta-1 < N)
			for(int j4=i+delta; j4<i+2*delta; ++j4)
				Right(i,1) += dst(Z(i,_), Z(j4,_), alpha);
	}

	//Update DLR
	for(int i = 1; i < minsize; ++i)
		for(int j = minsize; j < minsize+delta; ++j)
			DLR[minsize-1] = DLR[minsize-1] + dst(Z(i,_), Z(j,_), alpha);
	
	for(int s = minsize; s < N-delta; ++s){
		double r1 = Left(s,0); //Z[s] moved from the right to left segment so track affected distances
		double r2 = Right(s-delta,1); //Z[s-delta] has been removed from the sample under consideration
		double r3 = dst(Z(s,_),Z(s+delta,_),alpha); //Account for double counting in a1 and a2 below
		double a1 = Left(s+delta,1); //Z[s+delta] has been added to the sample under consideration
		double a2 = Right(s,0); //Z[s] moved from the right to left segment so track affected distances
		double a3 = dst(Z(s,_),Z(s-delta,_),alpha); //Account for double counting in r1 and r2 above
		DLR[s] = DLR[s-1] - r1 - r2 - r3 + a1 + a2 + a3;
	}

	//Calculaton of cumulative sum of distances for adjacent observations
	NumericVector cSum(N);
	std::fill(cSum.begin(), cSum.end(), 0.0);
	//cSum[i] = |Z[0]-Z[1]|^alpha + |Z[1]-Z[2]|^alpha + ... + |Z[i-1]-Z[i]|^alpha
	for(int i = 1;i < N; ++i)
		cSum[i] = cSum[i-1] + dst(Z(i,_), Z(i-1,_), alpha);

	if(verbose)
		Rcout<<"Pre-processing complete. Starting optimization."<<std::endl;

	//Solve for K=1
	std::vector<std::set<int> > testPoints(N);
	//testPoints[i] holds the pruned set of change point locations before i
	//Use std::set to make sure that there are no duplicate elements, plus have fast 
	//deletion of elements.
	// if K==1, only need to find best segmentation for entire series
	if(K==1){
	  int t = N-1;
	  std::set<int> cpSet; //Initially all locations are posible change point locations
	  for(int a=delta; a<=t-minsize; ++a)
	    cpSet.insert(a);
	  for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
	    //Iterate over possible change point locations for a given "ending position"
	    int s = *si;//change point location under consideration
	    int u = -1;//Only allowed to have 1 change point in this case (K=1)
	    cll = crr = delta*(delta-1)/2;//Initialize counters
	    clr = delta*delta;
	    
	    dll = DLL[s];//Initial within distnace sum for left segment
	    drr = DRR[s];//Initial within distance sum for right segment
	    dlr = DLR[s];//Initial between distance sum
	    //Augment distnace sums
	    dll += cSum[s-delta];
	    drr += ( cSum[t] - cSum[s+delta] );
	    cll += (s-delta-u);
	    crr += (t-s-delta);
	    
	    //Calculate test statistic
	    double stat = 2.0*dlr/(clr+0.0) - drr/(crr+0.0) - dll/(cll+0.0);
	    double num = (s-u)*(t-s)*(1.0);
	    double dnom = std::pow(t-u+0.0,2.0);
	    stat *= (num/dnom);
	    if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
	      FF(0,t) = stat;
	      A(0,t) = s;
	    }
	  }
	}else{
	  for(int t = 2*minsize-1; t < N; ++t){//Iterate over "ending position" for time series
	    std::set<int> cpSet; //Initially all locations are posible change point locations
	    for(int a=delta; a<=t-minsize; ++a)
	      cpSet.insert(a);
	    for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
	      //Iterate over possible change point locations for a given "ending position"
	      int s = *si;//change point location under consideration
	      int u = -1;//Only allowed to have 1 change point in this case (K=1)
	      cll = crr = delta*(delta-1)/2;//Initialize counters
	      clr = delta*delta;
	      
	      dll = DLL[s];//Initial within distnace sum for left segment
	      drr = DRR[s];//Initial within distance sum for right segment
	      dlr = DLR[s];//Initial between distance sum
	      //Augment distnace sums
	      dll += cSum[s-delta];
	      drr += ( cSum[t] - cSum[s+delta] );
	      cll += (s-delta-u);
	      crr += (t-s-delta);
	      
	      //Calculate test statistic
	      double stat = 2.0*dlr/(clr+0.0) - drr/(crr+0.0) - dll/(cll+0.0);
	      double num = (s-u)*(t-s)*(1.0);
	      double dnom = std::pow(t-u+0.0,2.0);
	      stat *= (num/dnom);
	      if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
	        FF(0,t) = stat;
	        A(0,t) = s;
	      }
	    }
	    //Update the set of possible change point locations for a given outer loop value
	    std::set<int> cpSetCopy = cpSet;
	    //Iterate over the copy and remove elements from the original
	    int a2 = t-minsize;		
	    int b2 = A(0,a2);
	    double V2;
	    if(b2 <= 0)
	      V2 = 0.0;
	    else
	      V2 = cSum[b2-1];
	    double cst2 = (a2-b2)*(t-a2)/std::pow(t-b2+0.0,2.0);
	    double stat2 = 2*DLR[a2]/(delta*delta);
	    double t1_2 = (DRR[a2]+cSum[t]-cSum[a2-delta])/(delta*(delta-1)/2+t-a2-delta);
	    double t2_2 = (DLL[a2]+cSum[a2-delta]-V2)/(delta*(delta-1)/2+a2-b2-delta);
	    stat2 -= (t1_2+t2_2);
	    stat2 *= cst2;
	    for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
	      int a = *si;
	      int b = A(0,a);
	      double V;
	      if(b <= 0)
	        V = 0.0;
	      else
	        V = cSum[b-1];
	      double cst = (a-b)*(t-a)/std::pow(t-b+0.0,2.0);
	      double stat = 2*DLR[a]/(delta*delta);
	      double t1 = (DRR[a]+cSum[t]-cSum[a-delta])/(delta*(delta-1)/2+t-a-delta);
	      double t2 = (DLL[a]+cSum[a-delta]-V)/(delta*(delta-1)/2+a-b-delta);
	      stat -= (t1+t2);
	      stat *= cst;
	      //Check pruning condition and remove elements as necessary
	      if(FF(0,a) + stat < FF(0,a2) + stat2)
	        cpSet.erase(a);
	    }
	    testPoints[t] = cpSet;
	  }
	}
	GOF[0] = FF(0,N-1);

	//If only 1 change point is to be found terminate early
	if(K == 1){
		if(verbose)
			Rcout<<"Finished optimization."<<std::endl;
		return wrap(List::create(_["number"]=1, _["estimates"]=A(0,N-1)+2, _["gofM"]=GOF));
	}

	//Solve for K > 1 cases
	//Since FF only has 2 rows we will use the variable 'flag' to help switch between the 0th and 1st row
	bool flag = true;

	for(int k=1; k<K; ++k){
		//Find change points
		// if (k+1)>=K, do not need to do pruning steps
		if((k+1)>=K){
		  int t = N-1;
		  FF(flag,t) = R_NegInf; //Reset goodness of fit value
		  std::set<int> cpSet = testPoints[t];
		  for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
		    int s = *si;
		    int u = A(k-1,s); //location of change point before s
		    cll = crr = delta*(delta-1)/2;//initialize counters
		    clr = delta*delta;
		    
		    dll = DLL[s]; //initial within distance sum for left segment
		    drr = DRR[s]; //initial within distnace sum for right segment
		    dlr = DLR[s]; //initial between distance sum
		    
		    //Augment distance sums
		    dll += cSum[s-delta];
		    if(u>0)
		      dll -= cSum[u-1];
		    drr += (cSum[t] - cSum[s+delta]);
		    cll += (s-delta-u);
		    crr += (t-s-delta);
		    
		    double stat = 2*dlr/(clr+0.0) - drr/(crr+0.0) - dll/(cll+0.0);
		    double num = (s-u)*(t-s)*1.0;
		    double dnom = std::pow(t-u+0.0, 2.0);
		    stat *= (num/dnom);
		    if(u>0) //Optimal partition cost if previous change points exists
		      stat += FF(1-flag,s);
		    if(stat > FF(flag,t)){
		      FF(flag,t) = stat;
		      A(k,t) = s;
		    }
		  }
		}else{
		  for(int t=2*minsize-1; t<N; ++t){
		    FF(flag,t) = R_NegInf; //Reset goodness of fit value
		    std::set<int> cpSet = testPoints[t];
		    for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
		      int s = *si;
		      int u = A(k-1,s); //location of change point before s
		      cll = crr = delta*(delta-1)/2;//initialize counters
		      clr = delta*delta;
		      
		      dll = DLL[s]; //initial within distance sum for left segment
		      drr = DRR[s]; //initial within distnace sum for right segment
		      dlr = DLR[s]; //initial between distance sum
		      
		      //Augment distance sums
		      dll += cSum[s-delta];
		      if(u>0)
		        dll -= cSum[u-1];
		      drr += (cSum[t] - cSum[s+delta]);
		      cll += (s-delta-u);
		      crr += (t-s-delta);
		      
		      double stat = 2*dlr/(clr+0.0) - drr/(crr+0.0) - dll/(cll+0.0);
		      double num = (s-u)*(t-s)*1.0;
		      double dnom = std::pow(t-u+0.0, 2.0);
		      stat *= (num/dnom);
		      if(u>0) //Optimal partition cost if previous change points exists
		        stat += FF(1-flag,s);
		      if(stat > FF(flag,t)){
		        FF(flag,t) = stat;
		        A(k,t) = s;
		      }
		    }
		    //Update the set of possible change point locations for a given set of outer loop values
		    std::set<int> cpSetCopy = cpSet;
		    int a2 = t-minsize;		
		    int b2 = A(k,a2);
		    double V2;
		    if(b2 <= 0)
		      V2 = 0.0;
		    else
		      V2 = cSum[b2-1];
		    double cst2 = (a2-b2)*(t-a2)*1.0/((t-b2)*(t-b2));
		    double stat2 = 2*DLR[a2]/(delta*delta);
		    double t1_2 = (DRR[a2]+cSum[t]-cSum[a2-delta])/(delta*(delta-1)/2+t-a2-delta);
		    double t2_2 = (DLL[a2]+cSum[a2-delta]-V2)/(delta*(delta-1)/2+a2-b2-delta);
		    stat2 -= (t1_2+t2_2);
		    stat2 *= cst2;
		    for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
		      int a = *si;
		      int b = A(k,a);
		      double V;
		      if(b <= 0)
		        V = 0;
		      else
		        V = cSum[b-1];
		      double cst = (a-b)*(t-a)*1.0/((t-b)*(t-b));
		      double stat = 2*DLR[a]/(delta*delta);
		      double t1 = (DRR[a]+cSum[t]-cSum[a-delta])/(delta*(delta-1)/2+t-a-delta);
		      double t2 = (DLL[a]+cSum[a-delta]-V)/(delta*(delta-1)/2+a-b-delta);
		      stat -= (t1+t2);
		      stat *= cst;
		      if(FF(flag,a) + stat < FF(flag,a2) + stat2)
		        cpSet.erase(a);
		    }
		    testPoints[t] = cpSet;
		  }
		}
		
		GOF[k] = FF(flag,N-1); //Obtain best goodness-of-fit value for the entire series
		flag = 1-flag;//Flip value of flag so that it now points to the alternate row
	}

	if(verbose)
		Rcout<<"Finished Optimization."<<std::endl;
	
	//Check stopping rule to determine the number of change points
	//according to point of inflection in GOF graph
	std::vector<double> POI;
	for(int i=1; i<(GOF.size()-1); ++i){
	  //regression for line 1
	  double x1_mean = 0.0;
	  for(int j=0; j<=i; ++j)
	    x1_mean += j;
	  x1_mean = x1_mean/(i+1);
	  double y1_mean = 0.0;
	  for(int j=0; j<=i; ++j)
	    y1_mean += GOF[j];
	  y1_mean = y1_mean/(i+1);
	  double beta1_num = 0.0, beta1_denom = 0.0;
	  for(int j=0; j<=i; ++j){
	    beta1_num += (j-x1_mean)*(GOF[j]-y1_mean);
	    beta1_denom += pow(j-x1_mean, 2);
	  }
	  double beta1 = beta1_num/beta1_denom;
	  double alpha1 = y1_mean - beta1*x1_mean;
	  //regression for line 2
	  double x2_mean = 0.0;
	  for(int j=i; j<GOF.size(); ++j)
	    x2_mean += j;
	  x2_mean = x2_mean/(GOF.size()-i);
	  double y2_mean = 0.0;
	  for(int j=i; j<GOF.size(); ++j)
	    y2_mean += GOF[j];
	  y2_mean = y2_mean/(GOF.size()-i);
	  double beta2_num = 0.0, beta2_denom = 0.0;
	  for(int j=i; j<GOF.size(); ++j){
	    beta2_num += (j-x2_mean)*(GOF[j]-y2_mean);
	    beta2_denom += pow(j-x2_mean, 2);
	  }
	  double beta2 = beta2_num/beta2_denom;
	  double alpha2 = y2_mean - beta2*x2_mean;
	  //fill in SSE
	  double SSE = 0.0;
	  for(int j=0; j<=i; ++j)
	    SSE += pow(alpha1 + beta1*j - GOF[j] ,2);
	  for(int j=i; j<GOF.size(); ++j)
	    SSE += pow(alpha2 + beta2*j - GOF[j] ,2);
	  POI.push_back(SSE);
	}
	//pick point of smallest SSE
	int k = -1;
	double SSE_inf = R_PosInf;
	for(int i=0; i<POI.size(); ++i){
	  if(POI[i]<SSE_inf){
	    SSE_inf = POI[i];
	    k = i;
	  }
	}
	
	//Find all optimal segmentations for differing numbers 
	//of change points.
	std::vector<std::vector<int> > cpLocs = find_locations(A);
	for (int i=0; i<K; ++i){
		transform(cpLocs[i].begin(), cpLocs[i].end(), cpLocs[i].begin(), std::bind(std::plus<int>(), placeholders::_1, 2));
	}
	std::vector<int> cps = cpLocs[k+1];

	return wrap(List::create(_["number"]=k+2, _["estimates"]=cps, _["gofM"]=GOF, _["cpLoc"]=cpLocs));

END_RCPP
}

SEXP eFastC( SEXP Z_, SEXP K_, SEXP minsize_, SEXP alpha_, SEXP verbose_){
BEGIN_RCPP
	
	//Convert SEXP variables to types Rcpp/C++ knows
	int K = as<int>(K_), minsize = as<int>(minsize_);
	double alpha = as<double>(alpha_);
	bool verbose = as<bool>(verbose_);
	NumericMatrix Z(Z_);
	
	int N = Z.nrow(); //Get number of observations
	int d = Z.ncol(); //Get number of dimensions

	// Reset K to be maximum number of possible segments
	if(N/minsize < K)
		K = N/minsize;

	//Since the update of the goodness-of-fit value with k change points only depends
	//upon the goodness-of-fit values with (k-1) change points we can use a 2xN matrix.
	//However, we will need a KxN matrix to store the change point locaitons.
	NumericMatrix FF(2,N), A(K,N); //Create matrices used to store goodness-of-fit values and 
	//change point locations

	//Goodness-of-fit values for when outer loop is at the end of the time series.
	//GOF[k] corresponds to the goodness-of-fit value  for the entire time series when
	//using k+1 change points.
	NumericVector GOF(K); 

	std::fill(FF.begin(), FF.end(), R_NegInf); //Fill FF matrix with negative infinity
	std::fill(A.begin(), A.end(), -1); //Fill A matrix with -1 because of C arrays being 0 indexed

	if(verbose)
		Rcout<<"Starting optimization."<<std::endl;

	//Solve for K=1
	std::vector<std::set<int> > testPoints(N);
	//testPoints[i] holds the pruned set of change point locations before i
	//Use std::set to make sure that there are no duplicate elements, plus have fast 
	//deletion of elements.
	// if K==1, only need to find best segmentation for entire series
	if(K==1){
	  if(verbose)
		Rcout<<"K = 1."<<std::endl;		
	  int t = N-1;
	  std::set<int> cpSet; //Initially all locations are posible change point locations
	  for(int a=minsize-1; a<=t-minsize; ++a)
	    cpSet.insert(a);
	  for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
	    //Iterate over possible change point locations for a given "ending position"
	    int s = *si;//change point location under consideration
	    int u = -1;//Only allowed to have 1 change point in this case (K=1)
	    double lr, ll, rr;
	    if(s==minsize-1){
	      NumericMatrix Z_us(s-u, d);
	      for (int i=u+1; i<=s; ++i)
	    	Z_us(i-u-1,_) = Z(i,_);
	      NumericMatrix Z_st(t-s, d);
	      for (int i=s+1; i<=t; ++i)
	    	Z_st(i-s-1,_) = Z(i,_);
	      lr = dist_XY(Z_us, Z_st, alpha);
	      ll = dist_X(Z_us, alpha);
	      rr = dist_X(Z_st, alpha);	    	
	    }else{
	      double lr_minus = 0.0;
	      double lr_plus = 0.0;
	      for (int i=0; i<s; i++){
	      	lr_minus += dst(Z(s,_), Z(i,_), alpha);
	      }
	      for (int i=s+1; i<=t; i++){
	      	lr_plus += dst(Z(s,_), Z(i,_), alpha);
	      }
	      lr = lr - lr_minus + lr_plus;
	      ll = ll + lr_minus;
	      rr = rr - lr_plus;
	    }
	    //Calculate test statistic
	    int n = s-u;
	    int m = t-s;
	    double stat = 2.0*lr/(n*m+0.0) - ll*2.0/(n*(n-1)+0.0) - rr*2.0/(m*(m-1)+0.0);
	    double num = (s-u)*(t-s)*(1.0);
	    double dnom = std::pow(t-u+0.0,2.0);
	    stat *= (num/dnom);
	    if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
	      FF(0,t) = stat;
	      A(0,t) = s;
	    }
	  }
	}else{
	  if(verbose)
		Rcout<<"K = 1."<<std::endl;			
	  for(int t = 2*minsize-1; t < N; ++t){//Iterate over "ending position" for time series
	    std::set<int> cpSet; //Initially all locations are posible change point locations
	    for(int a=minsize-1; a<=t-minsize; ++a)
	      cpSet.insert(a);
	    for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
	      //Iterate over possible change point locations for a given "ending position"
	      int s = *si;//change point location under consideration
	      int u = -1;//Only allowed to have 1 change point in this case (K=1)
	      double lr, ll, rr;
	      if(s==minsize-1){
	        NumericMatrix Z_us(s-u, d);
	        for (int i=u+1; i<=s; ++i)
	    	  Z_us(i-u-1,_) = Z(i,_);
	        NumericMatrix Z_st(t-s, d);
	        for (int i=s+1; i<=t; ++i)
	    	  Z_st(i-s-1,_) = Z(i,_);
	        lr = dist_XY(Z_us, Z_st, alpha);
	        ll = dist_X(Z_us, alpha);
	        rr = dist_X(Z_st, alpha);	    	
	      }else{
	        double lr_minus = 0.0;
	        double lr_plus = 0.0;
	        for (int i=0; i<s; i++){
	      	  lr_minus += dst(Z(s,_), Z(i,_), alpha);
	        }
	        for (int i=s+1; i<=t; i++){
	      	  lr_plus += dst(Z(s,_), Z(i,_), alpha);
	        }
	        lr = lr - lr_minus + lr_plus;
	        ll = ll + lr_minus;
	        rr = rr - lr_plus;
	      }
		  //Calculate test statistic
		  int n = s-u;
		  int m = t-s;
		  double stat = 2.0*lr/(n*m+0.0) - ll*2.0/(n*(n-1)+0.0) - rr*2.0/(m*(m-1)+0.0);
	      double num = (s-u)*(t-s)*(1.0);
	      double dnom = std::pow(t-u+0.0,2.0);
	      stat *= (num/dnom);
	      if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
	        FF(0,t) = stat;
	        A(0,t) = s;
	      }
	    }
	    //Update the set of possible change point locations for a given outer loop value
	    std::set<int> cpSetCopy = cpSet;
	    //Iterate over the copy and remove elements from the original
	    int a2 = t-minsize;		
	    int b2 = A(0,a2);
	    NumericMatrix Z_ba(a2-b2, d);
	    for (int i=b2+1; i<=a2; ++i)
	      Z_ba(i-b2-1,_) = Z(i,_);
	    NumericMatrix Z_at(t-a2, d);
	    for (int i=a2+1; i<=t; ++i)
	      Z_at(i-a2-1,_) = Z(i,_);
	    double lr = dist_XY(Z_ba, Z_at, alpha);
	    double ll = dist_X(Z_ba, alpha);
	    double rr = dist_X(Z_at, alpha);
	    //Calculate test statistic
	    int n = a2-b2;
	    int m = t-a2;
	    double stat2 = 2.0*lr/(n*m+0.0) - ll*2.0/(n*(n-1)+0.0) - rr*2.0/(m*(m-1)+0.0);
	    double cst2 = (a2-b2)*(t-a2)/std::pow(t-b2+0.0,2.0);
	    stat2 *= cst2;
	    for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
	      int a = *si;
	      int b = A(0,a);
	      NumericMatrix Z_ba(a-b, d);
	      for (int i=b+1; i<=a; ++i)
	    	Z_ba(i-b-1,_) = Z(i,_);
	      NumericMatrix Z_at(t-a, d);
	      for (int i=a+1; i<=t; ++i)
	    	Z_at(i-a-1,_) = Z(i,_);
	      double lr = dist_XY(Z_ba, Z_at, alpha);
	      double ll = dist_X(Z_ba, alpha);
	      double rr = dist_X(Z_at, alpha);
	      //Calculate test statistic
	      int n = a-b;
	      int m = t-a;
	      double stat = 2.0*lr/(n*m+0.0) - ll*2.0/(n*(n-1)+0.0) - rr*2.0/(m*(m-1)+0.0);
	      double cst = (a-b)*(t-a)/std::pow(t-b+0.0,2.0);
	      stat *= cst;
	      //Check pruning condition and remove elements as necessary
	      if(FF(0,a) + stat < FF(0,a2) + stat2)
	        cpSet.erase(a);
	    }
	    testPoints[t] = cpSet;
	  }
	}
	GOF[0] = FF(0,N-1);

	//If only 1 change point is to be found terminate early
	if(K == 1){
		if(verbose)
			Rcout<<"Finished optimization."<<std::endl;
		return wrap(List::create(_["number"]=1, _["estimates"]=A(0,N-1)+2, _["gofM"]=GOF));
	}

	//Solve for K > 1 cases
	//Since FF only has 2 rows we will use the variable 'flag' to help switch between the 0th and 1st row
	bool flag = true;

	for(int k=1; k<K; ++k){
	  	if(verbose)
		  Rcout<<"K = " << k+1 << '.'<<std::endl;				
		//Find change points
		// if (k+1)>=K, do not need to do pruning steps
		if((k+1)>=K){
		  int t = N-1;
		  FF(flag,t) = R_NegInf; //Reset goodness of fit value
		  std::set<int> cpSet = testPoints[t];
		  for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
		    int s = *si;
		    int u = A(k-1,s); //location of change point before s
	        NumericMatrix Z_us(s-u, d);
	        for (int i=u+1; i<=s; ++i)
	    	  Z_us(i-u-1,_) = Z(i,_);
	        NumericMatrix Z_st(t-s, d);
	        for (int i=s+1; i<=t; ++i)
	    	  Z_st(i-s-1,_) = Z(i,_);
	        double lr = dist_XY(Z_us, Z_st, alpha);
	        double ll = dist_X(Z_us, alpha);
	        double rr = dist_X(Z_st, alpha);
	        //Calculate test statistic
	        int n = s-u;
	        int m = t-s;
	        double stat = 2.0*lr/(n*m+0.0) - ll*2.0/(n*(n-1)+0.0) - rr*2.0/(m*(m-1)+0.0);
		    double num = (s-u)*(t-s)*1.0;
		    double dnom = std::pow(t-u+0.0, 2.0);
		    stat *= (num/dnom);
		    if(u>0) //Optimal partition cost if previous change points exists
		      stat += FF(1-flag,s);
		    if(stat > FF(flag,t)){
		      FF(flag,t) = stat;
		      A(k,t) = s;
		    }
		  }
		}else{
		  for(int t=2*minsize-1; t<N; ++t){		  	
		    FF(flag,t) = R_NegInf; //Reset goodness of fit value
		    std::set<int> cpSet = testPoints[t];
		    for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
		      int s = *si;
		      int u = A(k-1,s); //location of change point before s
	          NumericMatrix Z_us(s-u, d);
	          for (int i=u+1; i<=s; ++i){
	    	    Z_us(i-u-1,_) = Z(i,_);	          	
	    	  }
	          NumericMatrix Z_st(t-s, d);
	          for (int i=s+1; i<=t; ++i){
	    	    Z_st(i-s-1,_) = Z(i,_);
	          }
	          double lr = dist_XY(Z_us, Z_st, alpha);
	          double ll = dist_X(Z_us, alpha);
	          double rr = dist_X(Z_st, alpha);
	          //Calculate test statistic
	          int n = s-u;
	          int m = t-s;
	          double stat = 2.0*lr/(n*m+0.0) - ll*2.0/(n*(n-1)+0.0) - rr*2.0/(m*(m-1)+0.0);
		      double num = (s-u)*(t-s)*1.0;
		      double dnom = std::pow(t-u+0.0, 2.0);
		      stat *= (num/dnom);
		      if(u>0) //Optimal partition cost if previous change points exists
		        stat += FF(1-flag,s);
		      if(stat > FF(flag,t)){
		        FF(flag,t) = stat;
		        A(k,t) = s;
		      }
		    }
		    //Update the set of possible change point locations for a given set of outer loop values
		    std::set<int> cpSetCopy = cpSet;
		    int a2 = t-minsize;		
		    int b2 = A(k,a2);
	    	NumericMatrix Z_ba(a2-b2, d);
	    	for (int i=b2+1; i<=a2; ++i)
	      	  Z_ba(i-b2-1,_) = Z(i,_);
	    	NumericMatrix Z_at(t-a2, d);
	    	for (int i=a2+1; i<=t; ++i)
	      	  Z_at(i-a2-1,_) = Z(i,_);
	    	double lr = dist_XY(Z_ba, Z_at, alpha);
	    	double ll = dist_X(Z_ba, alpha);
	    	double rr = dist_X(Z_at, alpha);
	    	//Calculate test statistic
	    	int n = a2-b2;
	    	int m = t-a2;
	    	double stat2 = 2.0*lr/(n*m+0.0) - ll*2.0/(n*(n-1)+0.0) - rr*2.0/(m*(m-1)+0.0);
		    double cst2 = (a2-b2)*(t-a2)*1.0/((t-b2)*(t-b2));
		    stat2 *= cst2;
		    for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
		      int a = *si;
		      int b = A(k,a);
	    	  NumericMatrix Z_ba(a-b, d);
	    	  for (int i=b+1; i<=a; ++i)
	      	    Z_ba(i-b-1,_) = Z(i,_);
	    	  NumericMatrix Z_at(t-a, d);
	    	  for (int i=a+1; i<=t; ++i)
	      	    Z_at(i-a-1,_) = Z(i,_);
	    	  double lr = dist_XY(Z_ba, Z_at, alpha);
	    	  double ll = dist_X(Z_ba, alpha);
	    	  double rr = dist_X(Z_at, alpha);
	    	  //Calculate test statistic
	    	  int n = a-b;
	    	  int m = t-a;
	    	  double stat = 2.0*lr/(n*m+0.0) - ll*2.0/(n*(n-1)+0.0) - rr*2.0/(m*(m-1)+0.0);
		      double cst = (a-b)*(t-a)*1.0/((t-b)*(t-b));
		      stat *= cst;
		      if(FF(flag,a) + stat < FF(flag,a2) + stat2)
		        cpSet.erase(a);
		    }
		    testPoints[t] = cpSet;
		  }
		}
		
		GOF[k] = FF(flag,N-1); //Obtain best goodness-of-fit value for the entire series
		flag = 1-flag;//Flip value of flag so that it now points to the alternate row
	}

	if(verbose)
		Rcout<<"Finished Optimization."<<std::endl;
	
	//Check stopping rule to determine the number of change points
	//according to point of inflection in GOF graph
	std::vector<double> POI;
	for(int i=1; i<(GOF.size()-1); ++i){
	  //regression for line 1
	  double x1_mean = 0.0;
	  for(int j=0; j<=i; ++j)
	    x1_mean += j;
	  x1_mean = x1_mean/(i+1);
	  double y1_mean = 0.0;
	  for(int j=0; j<=i; ++j)
	    y1_mean += GOF[j];
	  y1_mean = y1_mean/(i+1);
	  double beta1_num = 0.0, beta1_denom = 0.0;
	  for(int j=0; j<=i; ++j){
	    beta1_num += (j-x1_mean)*(GOF[j]-y1_mean);
	    beta1_denom += pow(j-x1_mean, 2);
	  }
	  double beta1 = beta1_num/beta1_denom;
	  double alpha1 = y1_mean - beta1*x1_mean;
	  //regression for line 2
	  double x2_mean = 0.0;
	  for(int j=i; j<GOF.size(); ++j)
	    x2_mean += j;
	  x2_mean = x2_mean/(GOF.size()-i);
	  double y2_mean = 0.0;
	  for(int j=i; j<GOF.size(); ++j)
	    y2_mean += GOF[j];
	  y2_mean = y2_mean/(GOF.size()-i);
	  double beta2_num = 0.0, beta2_denom = 0.0;
	  for(int j=i; j<GOF.size(); ++j){
	    beta2_num += (j-x2_mean)*(GOF[j]-y2_mean);
	    beta2_denom += pow(j-x2_mean, 2);
	  }
	  double beta2 = beta2_num/beta2_denom;
	  double alpha2 = y2_mean - beta2*x2_mean;
	  //fill in SSE
	  double SSE = 0.0;
	  for(int j=0; j<=i; ++j)
	    SSE += pow(alpha1 + beta1*j - GOF[j] ,2);
	  for(int j=i; j<GOF.size(); ++j)
	    SSE += pow(alpha2 + beta2*j - GOF[j] ,2);
	  POI.push_back(SSE);
	}
	//pick point of smallest SSE
	int k = -1;
	double SSE_inf = R_PosInf;
	for(int i=0; i<POI.size(); ++i){
	  if(POI[i]<SSE_inf){
	    SSE_inf = POI[i];
	    k = i;
	  }
	}
	
	//Find all optimal segmentations for differing numbers 
	//of change points.
	std::vector<std::vector<int> > cpLocs = find_locations(A);
	for (int i=0; i<K; ++i){
		transform(cpLocs[i].begin(), cpLocs[i].end(), cpLocs[i].begin(), std::bind(std::plus<int>(), placeholders::_1, 2));
	}
	std::vector<int> cps = cpLocs[k+1];
	
	return wrap(List::create(_["number"]=k+2, _["estimates"]=cps, _["gofM"]=GOF, _["cpLoc"]=cpLocs));

END_RCPP
}

std::vector<double> MEAN_VAR(const std::vector<double>& X){
	//Calculate the mean and sample variane for observations in X
	//The variance is calculated in a single pass
	if(X.size() == 0){
		std::vector<double> ret;
		ret.push_back(0); ret.push_back(0);
		return ret;
	}

	double mean = *(X.begin()), var = 0;

	std::vector<double>::const_iterator i = X.begin();
	++i;
	int cnt = 2;
	for(; i!=X.end(); ++i, ++cnt){
		double dif = *i-mean;
		var = var*(cnt-2)/((double)cnt-1) + dif*dif/cnt;
		mean += dif/cnt;
	}
	
	std::vector<double> res;
	res.push_back(mean); res.push_back(var);
	return res;
}

double dst(const NumericVector& X, const NumericVector& Y, double alpha){
	//Calculate Euclidean distance between X and Y
	NumericVector res = X-Y;
	double ip = std::inner_product(res.begin(), res.end(), res.begin(),0.0);
	return std::pow( ip, alpha/2 );
}

double delta_sum(NumericMatrix& X, int a, int b, double alpha){
	//Determine the sum of the alpha distances between X[a], X[a+1], ..., X[b]
	double ret = 0.0;
	for(int i=a; i<b; ++i)
		for(int j=i+1; j<=b; ++j)
			ret += dst(X(i,_), X(j,_), alpha);
	return ret;
}

double dist_X(NumericMatrix& X, double alpha){
	int n = X.nrow();
	double ret = 0.0;
	for(int i=0; i<n-1; ++i)
		for(int j=i+1; j<n; ++j)
			ret += dst(X(i,_), X(j,_), alpha);
	return ret;	
}

double dist_XY(NumericMatrix& X, NumericMatrix& Y, double alpha){
	int n = X.nrow();
	int m = Y.nrow();
	double ret = 0.0;
	for(int i=0; i<n; ++i)
		for(int j=0; j<m; ++j)
			ret += dst(X(i,_), Y(j,_), alpha);
	return ret;	
}

std::vector<std::vector<int> > find_locations(const NumericMatrix& A){
	//Determine all of the optimal segmentations for 
	//differing numbers of change points
	std::vector<std::vector<int> > res;
	//Obtain dimensions for matrix A
	//N = number of observations, K = maximum number of fit change points
	int K = A.nrow(), N = A.ncol();
	//k+1 is the number of change points
	for(int k=0; k<K; ++k){
		int cp = A(k,N-1), k1 = k;
		std::vector<int> cps;
		do{
			cps.push_back(cp);
			--k1;
			if(k1 >= 0)
				cp = A(k1,cp);
		} while(k1 >= 0);
		sort(cps.begin(),cps.end());
		res.push_back(cps);
	}
	return res;
}

SEXP ksFastC( SEXP Z_, SEXP K_, SEXP minsize_, SEXP verbose_){
BEGIN_RCPP
	
	//Convert SEXP variables to types Rcpp/C++ knows
	int K = as<int>(K_), minsize = as<int>(minsize_);
	bool verbose = as<bool>(verbose_);
	NumericMatrix Z(Z_);
	
	int N = Z.nrow(); //Get number of observations

	// Reset K to be maximum number of possible segments
	if(N/minsize < K)
		K = N/minsize;

	//Since the update of the goodness-of-fit value with k change points only depends
	//upon the goodness-of-fit values with (k-1) change points we can use a 2xN matrix.
	//However, we will need a KxN matrix to store the change point locaitons.
	NumericMatrix FF(2,N), A(K,N); //Create matrices used to store goodness-of-fit values and 
	//change point locations

	//Goodness-of-fit values for when outer loop is at the end of the time series.
	//GOF[k] corresponds to the goodness-of-fit value  for the entire time series when
	//using k+1 change points.
	NumericVector GOF(K); 

	std::fill(FF.begin(), FF.end(), R_NegInf); //Fill FF matrix with negative infinity
	std::fill(A.begin(), A.end(), -1); //Fill A matrix with -1 because of C arrays being 0 indexed
	
	if(verbose)
		Rcout<<"Starting optimization."<<std::endl;
	
	//Solve for K=1
	std::vector<std::set<int> > testPoints(N);
	//testPoints[i] holds the pruned set of change point locations before i
	//Use std::set to make sure that there are no duplicate elements, plus have fast 
	//deletion of elements.
	std::vector<std::vector<pair<int,double> > > ksVal;
	ksVal.reserve(N-2*minsize+1);
	// ksVals contains the test statistics to be used for the next iteration of K
	// ksVals[t] holds the test statistics corresponding to change points in testPoints[t]
	// initializing
	std::set<int> cpSet;
	std::vector<double> vec_us_sort;
	vec_us_sort.reserve(N-minsize);
	std::vector<pair<double,int> > vec_st_sort;
	vec_st_sort.reserve(N-minsize);
	for(int t=2*minsize-1; t<N; ++t){//Iterate over "ending position" for time series
		if(t==2*minsize-1){// initialize at first t
			//Initially all locations are posible change point locations
			cpSet.insert(minsize-1);
			// initializing samples at first change point
			// s=minsize-1, u=-1
			for(int i=0; i<=minsize-1; ++i)
				vec_us_sort.push_back(Z[i]);
			for(int i=minsize; i<=t; ++i)
				vec_st_sort.push_back(std::make_pair(Z[i],i));
			// sorting. vec_st need to be resorted at each t
			std::sort(vec_us_sort.begin(),vec_us_sort.end());
			std::sort(vec_st_sort.begin(),vec_st_sort.end());
		}else{ // each iteration, add one element to series under consideration
			cpSet.insert(t-minsize);
			vec_st_sort.push_back(std::make_pair(Z[t],t));
			std::sort(vec_st_sort.begin(),vec_st_sort.end());
		}
		// initializing vector before and after change point
		std::vector<double> vec_us = vec_us_sort; 
		std::vector<pair<double,int> > vec_st_pair = vec_st_sort;
		std::vector<double> vec_st;
		vec_st.clear();
		vec_st.reserve(N-minsize);
		// initializing vector containing test stat for pruning
		std::vector<double> stat_prune;
		stat_prune.clear();
		stat_prune.reserve(cpSet.size());
		//Iterate over possible change point locations for a given "ending position" t
		for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
			int s = *si;//change point location under consideration
			if(s==minsize-1){// initialize at first change point
				// find the vector before and after change point
				for(int i=0; i<vec_st_sort.size(); ++i){
					vec_st.push_back(vec_st_sort[i].first);
				}
			}else{// each iteration, change point location will proceed by one, so vector before and after change point need to be updated
				// find index to be removed
				for(int to_erase=0; to_erase<vec_st_pair.size(); ++to_erase){
					if(vec_st_pair[to_erase].second==s){
						vec_us.push_back(vec_st_pair[to_erase].first);
						std::sort(vec_us.begin(),vec_us.end());
						vec_st_pair.erase(vec_st_pair.begin()+to_erase);
						vec_st.erase(vec_st.begin()+to_erase);
						break;
					}
				}
			}
			// if K==1, only need to calculate best segmentation for entire series
			if(K==1){
			  if(t==(N-1)){
			    //Calculate test statistic		
			    double stat = dist_ks(vec_us,vec_st);
			    if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
			      FF(0,t) = stat;
			      A(0,t) = s;
			    }
			  }
			}else{
			  //Calculate test statistic		
			  double stat = dist_ks(vec_us,vec_st);
			  if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
			    FF(0,t) = stat;
			    A(0,t) = s;
			  }
			  //Calculate test statistic for pruning
			  if(s<2*minsize-1){
			    stat_prune.push_back(stat);
			  }else{
			    int b = A(0,s);
			    vector<double> vec_bs(s-b);
			    for(int i=b+1; i<=s; ++i)
			      vec_bs[i-b-1] = Z[i];
			    std::sort(vec_bs.begin(),vec_bs.end());
			    stat_prune.push_back(dist_ks(vec_bs,vec_st));
			  }
			}
		}
		if(K>1){
  		//Update the set of possible change point locations for a given outer loop value
  		std::set<int> cpSetCopy = cpSet;
  		std::vector<pair<int,double> > ksVal_t; // element of ksVal
  		ksVal_t.clear();
  		ksVal_t.reserve(cpSet.size());
  		//Iterate over the original and remove elements from the copy
  		for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
  			int a = *si;
  			//Check pruning condition and remove elements as necessary
  			if(FF(0,a) + stat_prune[a-minsize+1] < FF(0,t-minsize) + stat_prune[t-2*minsize+1])
  				cpSetCopy.erase(a);
  			else
  				ksVal_t.push_back(std::make_pair(a,stat_prune[a-minsize+1]));
  		}
  		ksVal.push_back(ksVal_t); // may not contain change points minsize-1 and t-minsize
  		testPoints[t] = cpSetCopy;
		}
	}

	GOF[0] = FF(0,N-1);

	//If only 1 change point is to be found terminate early
	if(K == 1){
		if(verbose)
			Rcout<<"Finished optimization."<<std::endl;
		return wrap(List::create(_["number"]=1, _["estimates"]=A(0,N-1)+2, _["gofM"]=GOF));
	}

	//Solve for K > 1 cases
	//Since FF only has 2 rows we will use the variable 'flag' to help switch between the 0th and 1st row
	bool flag = true;
	
	for(int k=1; k<K; ++k){
		std::vector<std::vector<pair<int,double> > > ksValNew;
		ksValNew.clear();
		ksValNew.reserve(N-2*minsize+1);
		// if (k+1)>=K, do not need to do pruning steps
		if((k+1)>=K){
		  //Find change points
		  int t = N-1;
		  FF(flag,t) = R_NegInf; //Reset goodness of fit value
		  std::set<int> cpSet = testPoints[t];
		  // retrieve test stat calculated in (k-1)th iteration
		  std::vector<pair<int,double> > ksVal_t = ksVal[t-2*minsize+1];
		  // initializing vector containing test stat for pruning
		  std::vector<double> stat_prune;
		  stat_prune.clear();
		  stat_prune.reserve(cpSet.size());
		  for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
		    int s = *si;
		    int u = A(k-1,s); //location of change point before s
		    // retrieve test stat calculated in (k-1)th iteration
		    int j=0; // find location of saved test stat
		    double stat=0.0;				
		    if(ksVal_t.size()>0){
		      for(j=0; j<ksVal_t.size(); ++j){
		        if(ksVal_t[j].first==s)
		          break;
		      }
		    }		
		    //Calculate test statistic for pruning
		    stat = ksVal_t[j].second;
		    if(u>0){
		      //Optimal partition cost if previous change points exists
		      stat = stat + FF(1-flag,s);
		    } 
		    if(stat > FF(flag,t)){
		      FF(flag,t) = stat;
		      A(k,t) = s;
		    }
		  }
		}else{
		  //Trace all change points
		  std::set<int> locs_processing = testPoints[N-1];
		  std::set<int> locs;
		  locs.clear();
		  locs.insert(N-1);
		  for(int i=N-2; i>=2*minsize-1; --i){
		  	if(i==*locs_processing.rbegin()){
		  		locs.insert(i);
		  		locs_processing.erase(i);
		  		locs_processing.insert(testPoints[i].begin(), testPoints[i].end());
		  	}
		  }
		  //Find change points
		  for(int t=2*minsize-1; t<N; ++t){
		  	FF(flag,t) = R_NegInf; //Reset goodness of fit value
		  	if(t==*locs.begin()){
		  		locs.erase(t);
			    std::set<int> cpSet = testPoints[t];
			    // retrieve test stat calculated in (k-1)th iteration
			    std::vector<pair<int,double> > ksVal_t = ksVal[t-2*minsize+1];
			    // initializing vector containing test stat for pruning
			    std::vector<double> stat_prune;
			    stat_prune.clear();
			    stat_prune.reserve(cpSet.size());
			    for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
			      int s = *si;
			      int u = A(k-1,s); //location of change point before s
			      // retrieve test stat calculated in (k-1)th iteration
			      int j=0; // find location of saved test stat
			      double stat=0.0;				
			      if(ksVal_t.size()>0){
			        for(j=0; j<ksVal_t.size(); ++j){
			          if(ksVal_t[j].first==s)
			            break;
			        }
			      }		
			      //Calculate test statistic for pruning
		        stat = ksVal_t[j].second;
		        ksVal_t.erase(ksVal_t.begin()+j);
		        //Calculate test statistic for pruning
		        int b = A(k,s);
		        std::vector<double> vec_bs(s-b),vec_st(t-s);
		        for(int i=b+1; i<=s; ++i)
		          vec_bs[i-b-1] = Z[i];
		        for(int i=s+1; i<=t; ++i)
		          vec_st[i-s-1] = Z[i];
		        std::sort(vec_bs.begin(),vec_bs.end());
		        std::sort(vec_st.begin(),vec_st.end());
		        stat_prune.push_back(dist_ks(vec_bs,vec_st));
			      if(u>0) //Optimal partition cost if previous change points exists
			        stat += FF(1-flag,s);
			      if(stat > FF(flag,t)){
			        FF(flag,t) = stat;
			        A(k,t) = s;
			      }
			    }
		      //Update the set of possible change point locations for a given set of outer loop values
		      std::set<int> cpSetCopy = cpSet;
		      std::vector<pair<int,double> > ksValNew_t;
		      ksValNew_t.clear();
		      ksValNew_t.reserve(cpSet.size());
		      int prune_counter = 0;
		      for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
		        int a = *si;
		        if(FF(flag,a) + stat_prune[prune_counter] < FF(flag,t-minsize) + stat_prune.back()){
		          cpSet.erase(a);
		        }else{
		          ksValNew_t.push_back(std::make_pair(a,stat_prune[prune_counter]));
		        }
		        ++prune_counter;			
		      }
		      ksValNew.push_back(ksValNew_t); // may not contain change points minsize-1 and t-minsize
		      testPoints[t] = cpSet;
		  	}else{
		  		std::vector<pair<int,double> > ksValNew_t;
		  		ksValNew.push_back(ksValNew_t);
		  	}


		  }
		}
		
		ksVal.clear();
		ksVal = ksValNew;

		GOF[k] = FF(flag,N-1); //Obtain best goodness-of-fit value for the entire series
		flag = 1-flag;//Flip value of flag so that it now points to the alternate row
	}

	if(verbose)
		Rcout<<"Finished Optimization."<<std::endl;

	//Check stopping rule to determine the number of change points
	//according to point of inflection in GOF graph
	std::vector<double> POI;
	for(int i=1; i<(GOF.size()-1); ++i){
	  //regression for line 1
	  double x1_mean = 0.0;
	  for(int j=0; j<=i; ++j)
	    x1_mean += j;
	  x1_mean = x1_mean/(i+1);
	  double y1_mean = 0.0;
	  for(int j=0; j<=i; ++j)
	    y1_mean += GOF[j];
	  y1_mean = y1_mean/(i+1);
	  double beta1_num = 0.0, beta1_denom = 0.0;
	  for(int j=0; j<=i; ++j){
	    beta1_num += (j-x1_mean)*(GOF[j]-y1_mean);
	    beta1_denom += pow(j-x1_mean, 2);
	  }
	  double beta1 = beta1_num/beta1_denom;
	  double alpha1 = y1_mean - beta1*x1_mean;
	  //regression for line 2
	  double x2_mean = 0.0;
	  for(int j=i; j<GOF.size(); ++j)
	    x2_mean += j;
	  x2_mean = x2_mean/(GOF.size()-i);
	  double y2_mean = 0.0;
	  for(int j=i; j<GOF.size(); ++j)
	    y2_mean += GOF[j];
	  y2_mean = y2_mean/(GOF.size()-i);
	  double beta2_num = 0.0, beta2_denom = 0.0;
	  for(int j=i; j<GOF.size(); ++j){
	    beta2_num += (j-x2_mean)*(GOF[j]-y2_mean);
	    beta2_denom += pow(j-x2_mean, 2);
	  }
	  double beta2 = beta2_num/beta2_denom;
	  double alpha2 = y2_mean - beta2*x2_mean;
	  //fill in SSE
	  double SSE = 0.0;
	  for(int j=0; j<=i; ++j)
	    SSE += pow(alpha1 + beta1*j - GOF[j] ,2);
	  for(int j=i; j<GOF.size(); ++j)
	    SSE += pow(alpha2 + beta2*j - GOF[j] ,2);
	  POI.push_back(SSE);
	}
	//pick point of smallest SSE
	int k = -1;
	double SSE_inf = R_PosInf;
	for(int i=0; i<POI.size(); ++i){
	  if(POI[i]<SSE_inf){
	    SSE_inf = POI[i];
	    k = i;
	  }
	}
	
	//Find all optimal segmentations for differing numbers 
	//of change points.
	std::vector<std::vector<int> > cpLocs = find_locations(A);
	for (int i=0; i<K; ++i){
		transform(cpLocs[i].begin(), cpLocs[i].end(), cpLocs[i].begin(), std::bind(std::plus<int>(), placeholders::_1, 2));
	}
	std::vector<int> cps = cpLocs[k+1];
	
	return wrap(List::create(_["number"]=k+2, _["estimates"]=cps, _["gofM"]=GOF, _["cpLoc"]=cpLocs));

END_RCPP
}


SEXP ksFastC_delta( SEXP Z_, SEXP K_, SEXP minsize_, SEXP verbose_){
BEGIN_RCPP
	
	//Convert SEXP variables to types Rcpp/C++ knows
	int K = as<int>(K_), minsize = as<int>(minsize_);
	bool verbose = as<bool>(verbose_);
	NumericMatrix Z(Z_);
	
	int N = Z.nrow(); //Get number of observations

	// Reset K to be maximum number of possible segments
	if(N/minsize < K)
		K = N/minsize;

	//Since the update of the goodness-of-fit value with k change points only depends
	//upon the goodness-of-fit values with (k-1) change points we can use a 2xN matrix.
	//However, we will need a KxN matrix to store the change point locaitons.
	NumericMatrix FF(2,N), A(K,N); //Create matrices used to store goodness-of-fit values and 
	//change point locations

	//Goodness-of-fit values for when outer loop is at the end of the time series.
	//GOF[k] corresponds to the goodness-of-fit value  for the entire time series when
	//using k+1 change points.
	NumericVector GOF(K); 

	std::fill(FF.begin(), FF.end(), R_NegInf); //Fill FF matrix with negative infinity
	std::fill(A.begin(), A.end(), -1); //Fill A matrix with -1 because of C arrays being 0 indexed
	
	if(verbose)
		Rcout<<"Starting optimization."<<std::endl;
	
	//Solve for K=1
	std::vector<std::set<int> > testPoints(N);
	//testPoints[i] holds the pruned set of change point locations before i
	//Use std::set to make sure that there are no duplicate elements, plus have fast 
	//deletion of elements.
	std::vector<std::vector<pair<int,double> > > ksVal;
	ksVal.reserve(N-2*minsize+1);
	// ksVals contains the test statistics to be used for the next iteration of K
	// ksVals[t] holds the test statistics corresponding to change points in testPoints[t]
	// initializing
	std::vector<double> statvec(N);
	// statvec holds ks distance at each change point
	std::set<int> cpSet;
	for(int t=2*minsize-1; t<N; ++t){//Iterate over "ending position" for time series
		if(t==2*minsize-1){// initialize at first t
			//Initially all locations are posible change point locations
			cpSet.insert(minsize-1);
		}else{ // each iteration, add one element to series under consideration
			cpSet.insert(t-minsize);
		}
		// initializing vector containing test stat for pruning
		std::vector<double> stat_prune;
		stat_prune.clear();
		stat_prune.reserve(cpSet.size());
		//Iterate over possible change point locations for a given "ending position" t
		for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
			int s = *si;//change point location under consideration
			std::vector<double> vec_us(minsize);
			std::vector<double> vec_st(minsize);
			double stat;			
			if(t==2*minsize-1){
				for (int i=s-minsize+1; i<=s; ++i){
					vec_us[i-s+minsize-1] = Z[i];
				}
				for (int i=s+1; i<=t; ++i){
					vec_st[i-s-1] = Z[i];
				}
				std::sort(vec_us.begin(),vec_us.end());
				std::sort(vec_st.begin(),vec_st.end());				
			}
			// if K==1, only need to calculate best segmentation for entire series
			if(K==1){
			  if(t==(N-1)){
			    //Calculate test statistic
			  	if(s==t-minsize){
				  for (int i=s-minsize+1; i<=s; ++i){
					vec_us[i-s+minsize-1] = Z[i];
				  }
				  for (int i=s+1; i<=t; ++i){
					vec_st[i-s-1] = Z[i];
				  }
				  std::sort(vec_us.begin(),vec_us.end());
				  std::sort(vec_st.begin(),vec_st.end());			  		
			  	  stat = dist_ks(vec_us,vec_st);
			  	  statvec[s] = stat;	
			  	}else{
			  		stat = statvec[s];
			  	}
			    if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
			      FF(0,t) = stat;
			      A(0,t) = s;
			    }
			  }
			}else{
			  if(t==2*minsize-1){
			  	stat = dist_ks(vec_us,vec_st);
			  	statvec[s] = stat;	  	
			  }else{
			  	if(s==t-minsize){
				  for (int i=s-minsize+1; i<=s; ++i){
					vec_us[i-s+minsize-1] = Z[i];
				  }
				  for (int i=s+1; i<=t; ++i){
					vec_st[i-s-1] = Z[i];
				  }
				  std::sort(vec_us.begin(),vec_us.end());
				  std::sort(vec_st.begin(),vec_st.end());				  		
			  	  stat = dist_ks(vec_us,vec_st);
			  	  statvec[s] = stat;
			  	}else{
			  	  stat = statvec[s];
			  	}
			  }
			  //Calculate test statistic		
			  if(stat > FF(0,t)){//Update goodness-of-fit and location matrices
			    FF(0,t) = stat;
			    A(0,t) = s;
			  }
			  //Calculate test statistic for pruning
			  if(s<2*minsize-1){
			    stat_prune.push_back(stat);
			  }else{
			    stat_prune.push_back(statvec[s]);
			  }
			}
		}
		if(K>1){
  		//Update the set of possible change point locations for a given outer loop value
  		std::set<int> cpSetCopy = cpSet;
  		std::vector<pair<int,double> > ksVal_t; // element of ksVal
  		ksVal_t.clear();
  		ksVal_t.reserve(cpSet.size());
  		//Iterate over the original and remove elements from the copy
  		for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
  			int a = *si;
  			//Check pruning condition and remove elements as necessary
  			if(FF(0,a) + stat_prune[a-minsize+1] < FF(0,t-minsize) + stat_prune[t-2*minsize+1])
  				cpSetCopy.erase(a);
  			else
  				ksVal_t.push_back(std::make_pair(a,stat_prune[a-minsize+1]));
  		}
  		ksVal.push_back(ksVal_t); // may not contain change points minsize-1 and t-minsize
  		testPoints[t] = cpSetCopy;
		}
	}

	GOF[0] = FF(0,N-1);

	//If only 1 change point is to be found terminate early
	if(K == 1){
		if(verbose)
			Rcout<<"Finished optimization."<<std::endl;
		return wrap(List::create(_["number"]=1, _["estimates"]=A(0,N-1)+2, _["gofM"]=GOF));
	}

	//Solve for K > 1 cases
	//Since FF only has 2 rows we will use the variable 'flag' to help switch between the 0th and 1st row
	bool flag = true;
	
	for(int k=1; k<K; ++k){
		std::vector<std::vector<pair<int,double> > > ksValNew;
		ksValNew.clear();
		ksValNew.reserve(N-2*minsize+1);
		// if (k+1)>=K, do not need to do pruning steps
		if((k+1)>=K){
		  //Find change points
		  int t = N-1;
		  FF(flag,t) = R_NegInf; //Reset goodness of fit value
		  std::set<int> cpSet = testPoints[t];
		  // retrieve test stat calculated in (k-1)th iteration
		  std::vector<pair<int,double> > ksVal_t = ksVal[t-2*minsize+1];
		  // initializing vector containing test stat for pruning
		  std::vector<double> stat_prune;
		  stat_prune.clear();
		  stat_prune.reserve(cpSet.size());
		  for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
		    int s = *si;
		    int u = A(k-1,s); //location of change point before s
		    // retrieve test stat calculated in (k-1)th iteration
		    int j=0; // find location of saved test stat
		    double stat=0.0;				
		    if(ksVal_t.size()>0){
		      for(j=0; j<ksVal_t.size(); ++j){
		        if(ksVal_t[j].first==s)
		          break;
		      }
		    }		
		    //Calculate test statistic for pruning
		    stat = ksVal_t[j].second;
		    if(u>0){
		      //Optimal partition cost if previous change points exists
		      stat = stat + FF(1-flag,s);
		    } 
		    if(stat > FF(flag,t)){
		      FF(flag,t) = stat;
		      A(k,t) = s;
		    }
		  }
		}else{
		  //Trace all change points
		  std::set<int> locs_processing = testPoints[N-1];
		  std::set<int> locs;
		  locs.clear();
		  locs.insert(N-1);
		  for(int i=N-2; i>=2*minsize-1; --i){
		  	if(i==*locs_processing.rbegin()){
		  		locs.insert(i);
		  		locs_processing.erase(i);
		  		locs_processing.insert(testPoints[i].begin(), testPoints[i].end());
		  	}
		  }
		  //Find change points
		  for(int t=2*minsize-1; t<N; ++t){
		  	FF(flag,t) = R_NegInf; //Reset goodness of fit value
		  	if(t==*locs.begin()){
		  		locs.erase(t);
			    std::set<int> cpSet = testPoints[t];
			    // retrieve test stat calculated in (k-1)th iteration
			    std::vector<pair<int,double> > ksVal_t = ksVal[t-2*minsize+1];
			    // initializing vector containing test stat for pruning
			    std::vector<double> stat_prune;
			    stat_prune.clear();
			    stat_prune.reserve(cpSet.size());
			    for(std::set<int>::iterator si=cpSet.begin(); si!=cpSet.end(); ++si){
			      int s = *si;
			      int u = A(k-1,s); //location of change point before s
			      // retrieve test stat calculated in (k-1)th iteration
			      int j=0; // find location of saved test stat
			      double stat=0.0;				
			      if(ksVal_t.size()>0){
			        for(j=0; j<ksVal_t.size(); ++j){
			          if(ksVal_t[j].first==s)
			            break;
			        }
			      }		
			      //Calculate test statistic for pruning
		        stat = ksVal_t[j].second;
		        ksVal_t.erase(ksVal_t.begin()+j);
		        //Calculate test statistic for pruning
		        stat_prune.push_back(statvec[s]);
			      if(u>0) //Optimal partition cost if previous change points exists
			        stat += FF(1-flag,s);
			      if(stat > FF(flag,t)){
			        FF(flag,t) = stat;
			        A(k,t) = s;
			      }
			    }
		      //Update the set of possible change point locations for a given set of outer loop values
		      std::set<int> cpSetCopy = cpSet;
		      std::vector<pair<int,double> > ksValNew_t;
		      ksValNew_t.clear();
		      ksValNew_t.reserve(cpSet.size());
		      int prune_counter = 0;
		      for(std::set<int>::iterator si=cpSetCopy.begin(); si!=cpSetCopy.end(); ++si){
		        int a = *si;
		        if(FF(flag,a) + stat_prune[prune_counter] < FF(flag,t-minsize) + stat_prune.back()){
		          cpSet.erase(a);
		        }else{
		          ksValNew_t.push_back(std::make_pair(a,stat_prune[prune_counter]));
		        }
		        ++prune_counter;			
		      }
		      ksValNew.push_back(ksValNew_t); // may not contain change points minsize-1 and t-minsize
		      testPoints[t] = cpSet;
		  	}else{
		  		std::vector<pair<int,double> > ksValNew_t;
		  		ksValNew.push_back(ksValNew_t);
		  	}


		  }
		}
		
		ksVal.clear();
		ksVal = ksValNew;

		GOF[k] = FF(flag,N-1); //Obtain best goodness-of-fit value for the entire series
		flag = 1-flag;//Flip value of flag so that it now points to the alternate row
	}

	if(verbose)
		Rcout<<"Finished Optimization."<<std::endl;

	//Check stopping rule to determine the number of change points
	//according to point of inflection in GOF graph
	std::vector<double> POI;
	for(int i=1; i<(GOF.size()-1); ++i){
	  //regression for line 1
	  double x1_mean = 0.0;
	  for(int j=0; j<=i; ++j)
	    x1_mean += j;
	  x1_mean = x1_mean/(i+1);
	  double y1_mean = 0.0;
	  for(int j=0; j<=i; ++j)
	    y1_mean += GOF[j];
	  y1_mean = y1_mean/(i+1);
	  double beta1_num = 0.0, beta1_denom = 0.0;
	  for(int j=0; j<=i; ++j){
	    beta1_num += (j-x1_mean)*(GOF[j]-y1_mean);
	    beta1_denom += pow(j-x1_mean, 2);
	  }
	  double beta1 = beta1_num/beta1_denom;
	  double alpha1 = y1_mean - beta1*x1_mean;
	  //regression for line 2
	  double x2_mean = 0.0;
	  for(int j=i; j<GOF.size(); ++j)
	    x2_mean += j;
	  x2_mean = x2_mean/(GOF.size()-i);
	  double y2_mean = 0.0;
	  for(int j=i; j<GOF.size(); ++j)
	    y2_mean += GOF[j];
	  y2_mean = y2_mean/(GOF.size()-i);
	  double beta2_num = 0.0, beta2_denom = 0.0;
	  for(int j=i; j<GOF.size(); ++j){
	    beta2_num += (j-x2_mean)*(GOF[j]-y2_mean);
	    beta2_denom += pow(j-x2_mean, 2);
	  }
	  double beta2 = beta2_num/beta2_denom;
	  double alpha2 = y2_mean - beta2*x2_mean;
	  //fill in SSE
	  double SSE = 0.0;
	  for(int j=0; j<=i; ++j)
	    SSE += pow(alpha1 + beta1*j - GOF[j] ,2);
	  for(int j=i; j<GOF.size(); ++j)
	    SSE += pow(alpha2 + beta2*j - GOF[j] ,2);
	  POI.push_back(SSE);
	}
	//pick point of smallest SSE
	int k = -1;
	double SSE_inf = R_PosInf;
	for(int i=0; i<POI.size(); ++i){
	  if(POI[i]<SSE_inf){
	    SSE_inf = POI[i];
	    k = i;
	  }
	}
	
	//Find all optimal segmentations for differing numbers 
	//of change points.
	std::vector<std::vector<int> > cpLocs = find_locations(A);
	for (int i=0; i<K; ++i){
		transform(cpLocs[i].begin(), cpLocs[i].end(), cpLocs[i].begin(), std::bind(std::plus<int>(), placeholders::_1, 2));
	}
	std::vector<int> cps = cpLocs[k+1];
	
	return wrap(List::create(_["number"]=k+2, _["estimates"]=cps, _["gofM"]=GOF, _["cpLoc"]=cpLocs));

END_RCPP
}


double dist_ks(std::vector<double>& X, std::vector<double>& Y){
	// Sorted X and Y, get their sizes.
	int n = X.size();
	int m = Y.size();
	// Maximum absolute difference between empirical
	// cdf values.
	double max_diff = 0;
	// Store the current difference between empirical cdf values.
	double current_diff = 0;
	// Index of element in X or Y currently under consideration.
	int x_counter = 0;
	int y_counter = 0;
	// Determine the maximum difference betwen the empirical
	// cdf's of X and Y.
	// Iterate over X and Y comparing the elements at index
	// x_counter and y_counter.
	while (x_counter < n && y_counter < m){
		if (X[x_counter] < Y[y_counter]){
			current_diff = current_diff + 1.0/n;
			++x_counter;
		} else if (X[x_counter] > Y[y_counter]){
			current_diff = current_diff - 1.0/m;
			++y_counter;		
		} else {
			// If values are equal then we need to increment both
			// the x_counter and y_counter.
			double s = X[x_counter];
			while (x_counter < n && X[x_counter] == s){
				current_diff = current_diff + 1.0/n;
				++x_counter;
			}
			while (y_counter < m && Y[y_counter] == s){
				current_diff = current_diff - 1.0/m;
				++y_counter;
			}
		}
		// Update the current value of max_diff.
		max_diff = (fabs(current_diff) > max_diff) ?
		fabs(current_diff) : max_diff;
	}
	// Since (x_counter < n && y_counter < m) is false then either
	// x_counter == n or y_counter == m. This while loop deals with the
	// first case. The following with the later.
	while (y_counter < m) {
		current_diff = current_diff - 1.0/m;
		++y_counter;
		max_diff = (fabs(current_diff) > max_diff) ?
		fabs(current_diff) : max_diff;
	}
	while (x_counter < n) {
		current_diff = current_diff + 1.0/n;
		++x_counter;
		max_diff = (fabs(current_diff) > max_diff) ?
		fabs(current_diff) : max_diff;
	}
	return max_diff*n*m/std::pow(n+m,2);
}
