
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP getWithin( SEXP alpha_, SEXP X_) ;
SEXP getBetween( SEXP alpha_, SEXP X_, SEXP Y_) ;
SEXP splitPointC( SEXP s_, SEXP e_, SEXP D_, SEXP min_size_) ;
SEXP getBounds(SEXP n_, SEXP lvl_, SEXP eps_) ;
}

// definition

SEXP getWithin( SEXP alpha_, SEXP X_ ){
BEGIN_RCPP
try{
NumericMatrix X(X_);//matrix of points used
		double alpha=as<double>(alpha_);
		double ret = 0.0;
		int n = X.nrow();
		for(int i=0;i<n;++i)
			for(int j=0;j<n;++j)
				ret += std::pow((sqrt(sum((X(i,_)-X(j,_))*(X(i,_)-X(j,_))))), alpha);
		return wrap(ret/(n*n));		
	}
catch(std::exception& ex){
		forward_exception_to_r(ex);
	}
	catch(...){
		Rf_error("unknown C++ exception");
	}
END_RCPP
}


SEXP getBetween( SEXP alpha_, SEXP X_, SEXP Y_ ){
BEGIN_RCPP

try{
		NumericMatrix X(X_), Y(Y_);
		double alpha=as<double>(alpha_);
		double ret = 0.0;
		int n = X.nrow(), m = Y.nrow();
		for(int i=0;i<n;++i)
			for(int j=0;j<m;++j)
				ret += std::pow((sqrt(sum((X(i,_)-Y(j,_))*(X(i,_)-Y(j,_))))), alpha);				
		return wrap(2*ret/(n*m));	
	}
	catch(std::exception& ex){
		forward_exception_to_r(ex);
	}
	catch(...){
		Rf_error("unknown C++ exception");
	}

END_RCPP
}


SEXP splitPointC( SEXP s_, SEXP e_, SEXP D_, SEXP min_size_ ){
BEGIN_RCPP

//This impementation takes at most O(n^2) time to find each change point.
//Thus if k change points are found, our algorithm is O(kn^2).
using namespace Rcpp;
//used to sum an individual row/colum of a matrix
#define SUM(A) std::accumulate(A.begin() , A.end() , 0.0)

NumericVector best = NumericVector::create(-1.0,R_NegInf);
int e = as<int>(e_), s = as<int>(s_), min_size = as<int>(min_size_);//ending index, starting index, minimum size

NumericMatrix D(D_);//the distance matrix
e = e-s+1;//now represents the number of data points

double t1=min_size, t2=min_size<<1;//tau1 and tau2
NumericMatrix cut1 = D(Range(0,t1-1),Range(0,t1-1));
NumericMatrix cut2 = D(Range(t1,t2-1),Range(t1,t2-1));
NumericMatrix cut3 = D(Range(0,t1-1),Range(t1,t2-1));

double A = SUM(cut1)/2, B1=SUM(cut2)/2, AB1=SUM(cut3);
double tmp= 2*AB1/((t2-t1)*(t1)) - 2*B1/((t2-t1-1)*(t2-t1)) - 2*A/((t1-1)*(t1));
tmp *= (t1*(t2-t1)/t2);
if(tmp > best[1]){//update if better location is found
	best[0] = t1+s;
	best[1] = tmp;
}

t2+=1;

NumericVector B(e+1,B1), AB(e+1,AB1);//might be wasting a little space because of
//elements at the beginning that are not used
for(;t2<=e;++t2){//update between and within (for the right sample) distances
	B[t2] = B[t2-1] + SUM(D(Range(t2-1,t2-1),Range(t1,t2-2)));
	AB[t2] = AB[t2-1] + SUM(D(Range(t2-1,t2-1),Range(0,t1-1)));
	tmp = 2*AB[t2]/((t2-t1)*(t1))-2*B[t2]/((t2-t1-1)*(t2-t1))-2*A/((t1)*(t1-1));
	tmp *= (t1*(t2-t1)/t2);
	if(tmp > best[1]){//update if better location is found
		best[0] = t1+s;
		best[1] = tmp;
	}
}

t1+=1;

for(;;++t1){//iterate over possible change point locations (t1)
	t2=t1+min_size;//skip tau2 to its smallest allowed value
	if(t2>e)//remaining interval is too small to fit another cluster
		break;
	double addA = SUM(D(Range(t1-1,t1-1),Range(0,t1-2)));
	A+=addA;//update within distance for left cluster
	double addB = SUM(D(Range(t1-1,t1-1),Range(t1,t2-2)));
	for(;t2<=e;++t2){//iterate over possible ending locations for right cluster (t2)
		addB += D(t1-1,t2-1);
		B[t2]-=addB;//update within disance for right cluster
		AB[t2]+=(addB-addA);//update between cluster distance
		tmp = 2*AB[t2]/((t2-t1)*(t1))-2*B[t2]/((t2-t1-1)*(t2-t1))-2*A/((t1-1)*(t1));
		tmp *= (t1*(t2-t1)/t2);//new test statistic
		if(tmp > best[1]){
			best[0] = t1+s;
			best[1] = tmp;
		}
	}
}
return wrap(best);

END_RCPP
}



SEXP getBounds(SEXP n_, SEXP lvl_, SEXP eps_){
BEGIN_RCPP
	//n_ = number of boundary values to generate
	//lvl_ = desired significance level
	//eps_ = epsilon spending vector
	int n = as<int>(n_), start = 0, Uat = 2, Lat = -1;
	double alpha = as<double>(lvl_), errU = 0.0, errL = 0.0;
	std::vector<int> U(n+1,0), L(n+1,0);//vectors to hold upper and lower bounds
	std::vector<double> prob(n+1,0), eps = as<std::vector<double> >(eps_);//prob = probability of hitting a boundary
	prob[0] = 1-alpha; prob[1] = alpha;
	U[0] = 2; L[0] = -1;//impossible to hit boundary on first attempt
	std::vector<double>::iterator vi,vi2;
	for(int i = 0; i<n; ++i){//sequentially calculate boundary values
		prob[Uat] = prob[Uat-1]*alpha;
		for(vi=prob.begin()+Uat-1,vi2=vi-1; vi!=prob.begin();--vi,--vi2)
			(*vi) = (*vi)*(1-alpha)+(*vi2)*alpha;
		prob[Lat+1] *= (1-alpha);
		while(prob[Uat]+errU <= eps[i+1]){//while still have spending room decrease upper bound value
			errU += prob[Uat];
			--Uat;
		}
		while(prob[Lat]+errL <= eps[i+1]){//while still have spending room increase lower bound value
			errL += prob[Uat];
			++Lat;
		}
		++Uat;
		L[i+1] = Lat+start;
		U[i+1] = Uat+start;
		if(Lat >= 0)//reset to -1, can reuse calculations and above code
			for(vi=prob.begin()+Lat+1,vi2=prob.begin();vi!=prob.begin()+Uat && vi!=prob.end();++vi,++vi2)
				(*vi2) = (*vi);
			start += Lat+1;
			Uat -= Lat+1;
			Lat = -1;
	}
	return List::create(_["u"]=U, _["l"]=L);
END_RCPP
}
