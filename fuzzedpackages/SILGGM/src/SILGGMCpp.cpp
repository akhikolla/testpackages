#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace std;
using namespace Rcpp;


/*
*** By Rong, 2017.10.09

*** Main references:
1. Eddelbuettel, D. et al. (2011) Rcpp: Seamless R and C++ integration. Journal of Statistical Software, 40, 1-18.
2. Friedman, J. et al. (2008) Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9, 432-441.
3. Friedman, J. et al. (2010) Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33, 1-22.
4. Jankova, J. and van de Geer, S. (2015) Confidence intervals for high-dimensional inverse covariance estimation. Electronic Journal of Statistics, 9, 1205-1229.
5. Jankova, J. and van de Geer, S. (2017) Honest confidence regions and optimality in high-dimensional precision matrix estimation. Test, 26, 143-162.
6. Liu, W. (2013) Gaussian graphical model estimation with false discovery rate control. The Annals of Statistics, 41, 2948-2978.
7. Ren, Z. et al. (2015) Asymptotic normality and optimalities in estimation of large Gaussian graphical models. The Annals of Statistics, 43, 991-1026.
8. Shannon, P. et al. (2003) Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Research, 13, 2498-2504.
9. Wang, T. et al. (2016) FastGGM: an efficient algorithm for the inference of gaussi-an graphical model in biological networks. PLoS Computational Biology, 12, e1004755.
10. Witten, D. M. et al. (2011) New insights and faster computations for the graphical lasso. Journal of Computational and Graphical Statistics, 20, 892-900.


*** Function: Statistical inference of large-scale Gaussian graphical model in gene networks


*/

/* ****** Sub-function, fast Lasso regression ******* */
NumericVector FastLasso_Rcpp(NumericVector ipy, NumericMatrix ipx, double lambda, size_t N){
  double tol = 1e-3; //threshold for stopping
  size_t p = ipy.size();
  double dbm;  //maximum of beta difference for while loop
  NumericVector beta(p); //initialize output beta vector with length p and filled with 0.0
  NumericVector gc(p); //initialize grandient components

  //update beta vector with coordinate descent, covariance updates
  do{
    dbm = 0;
    for(size_t j = 0; j < p; j++){
      double z = (ipy[j] - gc[j])/N + beta[j];
      double beta_tmp = max(0.0, z - lambda) - max(0.0, -z - lambda);
      double diffBeta = beta_tmp - beta[j];
      double diffabs = abs(diffBeta);
      if (diffabs > 0){
        beta[j] = beta_tmp;
        for (size_t i = 0; i < p; i++){
          gc[i] = gc[i] + ipx(i,j) * diffBeta;
        }
        dbm = max(dbm, diffabs);
      }
    }
  }
  while(dbm >= tol);

  return beta;
}

/* ****** Sub-function, fast Lasso regression used for warm start ******* */
List FastLasso_Rcpp_v2(NumericVector beta, NumericVector gc, NumericVector ipy, NumericMatrix ipx, double lambda, size_t N){
  double tol = 1e-3; //threshold for stopping
  size_t p = ipy.size();
  double dbm;  //maximum of beta difference for while loop

  //update beta vector with coordinate descent, covariance updates
  do{
    dbm = 0;
    for(size_t j = 0; j < p; j++){
      double z = (ipy[j] - gc[j])/N + beta[j];
      double beta_tmp = max(0.0, z - lambda) - max(0.0, -z - lambda);
      double diffBeta = beta_tmp - beta[j];
      double diffabs = abs(diffBeta);
      if (diffabs > 0){
        beta[j] = beta_tmp;
        for (size_t i = 0; i < p; i++){
          gc[i] = gc[i] + ipx(i,j) * diffBeta;
        }
        dbm = max(dbm, diffabs);
      }
    }
  }
  while(dbm >= tol);

  return List::create(_["beta"] = beta, _["gc"] = gc);
}



/* ***** return index of minimum element ***** */
double min_position(NumericVector x){
	 NumericVector::iterator it = min_element(x.begin(), x.end());
	 double position = it - x.begin();
	 return position;
}

/* ***** sort the vector in descending order ***** */
NumericVector sortDescending(NumericVector x){
	sort(x.begin(), x.end(), greater<double>());
	return x;
}



/* ***** Function 1 ******** */
// [[Rcpp::export]]
List SILGGMCpp(NumericMatrix x, Nullable<CharacterVector> method=R_NilValue, Nullable<double> lambda=R_NilValue, bool global = false, Nullable<NumericVector> alpha=R_NilValue, Nullable<double> ndelta=R_NilValue, Nullable<NumericMatrix> true_graph=R_NilValue){
  // Stop threshold and parameter for edge decision
  double tol = 1e-3;

  size_t n = x.nrow();
  size_t p = x.ncol();

  CharacterVector rcpp_method = CharacterVector::create();
  if(method.isNull()){
	  rcpp_method = "D-S_NW_SL";
	  Rcout << "Use default method 'D-S_NW_SL'" << endl;
  }
  else{
	  rcpp_method = as<CharacterVector>(method);
	  Rcout << "Use method '" << rcpp_method << "'" << endl;
  }

  if(rcpp_method[0] != "D-S_NW_SL" && rcpp_method[0] != "D-S_GL" && rcpp_method[0] != "GFC_SL" && rcpp_method[0] != "GFC_L" && rcpp_method[0] != "B_NW_SL"){
	  stop("The method is not available, please try a correct method: 'D-S_NW_SL', 'D-S_GL', 'B_NW_SL', 'GFC_SL' or 'GFC_L'");
  }

  NumericVector rcpp_level = NumericVector::create();
  if(alpha.isNull()){
  	rcpp_level = NumericVector::create(0.05, 0.1);
  }
  else{
  	rcpp_level = as<NumericVector>(alpha);
  }

  size_t alpha_size = rcpp_level.size();


  if(rcpp_method[0] == "D-S_NW_SL"){
  NumericMatrix Omega(p,p);
  NumericMatrix partialCor(p,p);
  NumericMatrix CI_low_Omega(p,p);
  NumericMatrix CI_high_Omega(p,p);
  NumericMatrix z_score(p,p);
  NumericMatrix p_Omega(p,p);
  NumericVector FDR(alpha_size);
  NumericVector Power(alpha_size);
  NumericVector threshold(alpha_size);
  List global_decision(alpha_size);

  //penalty parameter for L1 norm of betas
  double rcpp_lambda;
  if(lambda.isNull()){
    rcpp_lambda = sqrt(2 * log(p/sqrt(n))/n);
    Rcout << "Use default lambda = sqrt(2*log(p/sqrt(n))/n)" << endl;
    Rcout << "In this case, lambda = " << rcpp_lambda << endl;
  }
  else{
	rcpp_lambda = as<double>(lambda);
	Rcout << "In this case, lambda = " << rcpp_lambda << endl;
  }

  // Center matrix
  Rcout << "Center each column." << endl;
  NumericMatrix x_cent(n,p);
  for(size_t i = 0; i < p; i++){
    x_cent(_,i) = x(_,i) - mean(x(_,i));
  }

  // Standardize matrix
  Rcout << "Standardize each column." << endl;
  NumericVector x_l2norm(p);
  for(size_t i = 0; i < p; i++){
    x_l2norm[i] = sqrt(sum(pow(x_cent(_,i),2)));
    if(x_l2norm[i] == 0){
      stop("some variables have 0 variance, please remove them at first.");
   }
  }

  NumericMatrix x_stan(n,p);
  for(size_t i = 0; i < p; i++){
    x_stan(_,i) = x_cent(_,i) / x_l2norm[i] * sqrt(n);
  }


  // Pre-calculate inner product matrixes
  Rcout << "Pre-calculate inner product matrixes." << endl;
  NumericMatrix IP_YX(p,p);
  for(size_t i = 0; i < p; i++){
    for(size_t j = 0; j < p; j++){
      IP_YX(i,j) = inner_product(x_cent(_,i).begin(), x_cent(_,i).end(), x_stan(_,j).begin(), 0.0);
    }
  }
  NumericMatrix IP_XX(p,p);
  NumericMatrix IP_YX_origin(p,p);
  for(size_t i = 0; i < p; i++){
    IP_XX(i,_) = IP_YX(i,_) / x_l2norm[i] * sqrt(n);
    IP_YX_origin(_,i) = IP_YX(_,i) * x_l2norm[i] / sqrt(n);
  }


  // Calculate scaled Lasso for each variable
  Rcout << "Calculate scaled Lasso for each variable." << endl;
  NumericMatrix beta_initial(p,p);
  NumericMatrix epsi_pre(n,p);
  NumericVector homega(p);
  NumericMatrix inner(p,p);
  NumericMatrix inner1(p,p);
  NumericMatrix precision(p,p);
  NumericMatrix precision1(p,p);

  for(size_t i = 0; i < p; i++){
    if(i % 100 == 0){
      Rcout <<"scaled Lasso for variable " << (i + 1) << endl;
    }
    double sigma = 1.0;
    NumericVector beta(p-1);
    NumericVector epsi(n);
    double reg, diffBeta, sigma_tmp, diffSigma;

    // Extract IP_XX[-i,-i] and IP_YX[i,-i]
    NumericMatrix tmpxx(p-1,p-1);
    NumericVector tmpyx(p-1);
    size_t t = 0;
    for (size_t k = 0; k < p; k++){
      if (k != i){
        tmpyx[t] = IP_YX(i,k);
        size_t t1 = 0;
        for (size_t l = 0; l < p; l++){ // l is index of row, k is index of col
          if (l != i){
            tmpxx(t1,t) = IP_XX(l,k);
            t1++;
          }
        }
        t++;
      }
    }

    // Extract x_stan[,-i]
    NumericMatrix tmp_x(n, p-1);
    t=0;
    for (size_t k = 0; k < p; k++){
      if (k != i){
        tmp_x(_,t) = x_stan(_,k);
        t++;
      }
    }

    // Extract x_l2norm[,-i]
    NumericVector x_l2norm_new(p-1);
    t=0;
    for (size_t k = 0; k < p; k++){
      if (k != i){
    	x_l2norm_new[t] = x_l2norm[k];
    	t++;
      }
    }


    //Scaled_lasso
    size_t iter = 0;
    NumericVector tmp_y(n);
    do {
      //Update beta when fix sigma
      reg = sigma * rcpp_lambda;
      NumericVector beta_tmp = FastLasso_Rcpp(tmpyx, tmpxx, reg, n);
      diffBeta = max(abs(beta_tmp - beta));
      beta = beta_tmp;
      //Update sigma when fix beta

      size_t non_zero_count = 0;
      for (size_t k = 0; k < p-1; k++){
    	  if(beta[k] != 0){
    		  non_zero_count++;
    	  }
      }

      NumericMatrix tmp_x_nonzero(n,non_zero_count);
      NumericVector beta_nonzero(non_zero_count);

      size_t t = 0;
      for (size_t k = 0; k < p-1; k++){
    	  if(beta[k] != 0){
    		  tmp_x_nonzero(_,t) = tmp_x(_,k);
    		  beta_nonzero[t] = beta[k];
    		  t++;
    	  }
      }

      for (size_t k = 0; k < n; k++){
        tmp_y[k] = inner_product(tmp_x_nonzero(k,_).begin(), tmp_x_nonzero(k,_).end(), beta_nonzero.begin(), 0.0);
      } //Multiplication of x.stan[,-i] and beta
      epsi = x_cent(_,i) - tmp_y;
      sigma_tmp = sqrt(sum(pow(epsi,2))/n);
      diffSigma = abs(sigma_tmp - sigma);
      sigma = sigma_tmp;
      iter++;
    }
    while((diffSigma >= tol || diffBeta >= tol) && iter < 10); //Stop iteration


    for (size_t k = 0; k < p; k++){
    	if(i < k){
    		beta_initial(i,k) = beta[k-1]/x_l2norm_new[k-1]*sqrt(n);
    	}
    	if(i > k){
    		beta_initial(i,k) = beta[k]/x_l2norm_new[k]*sqrt(n);
    	}
    }

    homega[i]= 1/(pow(sigma,2)+rcpp_lambda*sigma*sum(abs(beta_initial(i,_))));
    precision(i,i) = homega[i];
    precision1(i,i) = homega[i];

  }

    //Bias Correction

    for (size_t i = 0; i < p-1; i++){
    	for(size_t j=i+1; j < p; j++){

            precision(i,j) = -beta_initial(i,j) * homega[i];
            precision(j,i) = -beta_initial(j,i) * homega[j];
            precision1(i,j) = precision(j,i);
            precision1(j,i) = precision(i,j);

   	  }
    }


    for(size_t j = 0; j < p; j++){
  	  size_t non_zero_count = 0;
  	  NumericVector beta = precision1(_,j);
  	  for (size_t k = 0; k < p; k++){
  			 if(beta[k] != 0){
  			     non_zero_count++;
  			 }
  	  }
  	  NumericMatrix S_nonzero(p,non_zero_count);
  	  NumericVector beta_nonzero(non_zero_count);

  	  size_t t = 0;
  	  for (size_t k = 0; k < p; k++){
  		    if(beta[k] != 0){
  		    	S_nonzero(_,t) = IP_YX_origin(_,k)/n;
  		    	beta_nonzero[t] = beta[k];
  		    	t++;
  		    }
  	  }

  	  for(size_t i = 0; i < p; i++){
  	     inner(i,j) = inner_product(S_nonzero(i,_).begin(), S_nonzero(i,_).end(), beta_nonzero.begin(), 0.0);
  	  }
    }

    for(size_t i = 0; i < p; i++){
  	  size_t non_zero_count = 0;
  	  NumericVector beta = precision(i,_);
  	  for (size_t k = 0; k < p; k++){
  		  if(beta[k] != 0){
  			  non_zero_count++;
  		  }
  	  }
  	  NumericMatrix inner_nonzero(non_zero_count, p);
  	  NumericVector beta_nonzero(non_zero_count);

  	  size_t t = 0;
  		  for (size_t k = 0; k < p; k++){
  			    if(beta[k] != 0){
  			    	inner_nonzero(t,_) = inner(k,_);
  			    	beta_nonzero[t] = beta[k];
  			    	t++;
  			    }
  		  }


  	  for(size_t j = 0; j < p; j++){
  		  inner1(i,j) = inner_product(beta_nonzero.begin(), beta_nonzero.end(), inner_nonzero(_,j).begin(),0.0);
  		  Omega(i,j) = precision1(i,j)+precision(i,j) - inner1(i,j);
  	  }
    }

  for(size_t i = 0; i < p-1; i++){
	  for(size_t j =i+1; j < p; j++){
		  partialCor(i,j) = -Omega(i,j) / sqrt(Omega(i,i)*Omega(j,j));
		  double std_new = sqrt((precision1(i,i)*precision1(j,j)+pow((precision1(i,j)+precision1(j,i))/2,2))/n);
		  z_score(i,j) = Omega(i,j)/std_new;
		  double z_95CI = 1.96;
		  CI_low_Omega(i,j) = Omega(i,j) - z_95CI * std_new;
		  CI_high_Omega(i,j) = Omega(i,j) + z_95CI * std_new;
		  double abs_zscore = abs(z_score(i,j));
		  p_Omega(i,j) = 2 * Rf_pnorm5(-abs_zscore, 0.0, 1.0, 1, 0);

		  partialCor(j,i) = partialCor(i,j);
		  z_score(j,i) = z_score(i,j);
		  CI_low_Omega(j,i) = CI_low_Omega(i,j);
		  CI_high_Omega(j,i) = CI_high_Omega(i,j);
		  p_Omega(j,i) = p_Omega(i,j);
	  }
  }

  if(global == true){

	    Rcout <<"Perform global inference." << endl;
	  	Rcout << "Use pre-specified level(s): " << rcpp_level << endl;

	    //Take out the upper triangular of test statistic
	    size_t t = 0;
	    NumericVector upper_stat(p*(p-1)/2);
	    for (size_t i = 1; i < p; i++){
	    	for (size_t j = 0; j < i; j++){
	    		upper_stat[t]=z_score(j,i);
	    		t++;
	    	}
	    }

	   //Take out the upper triangular of true graph if available
	    NumericMatrix rcpp_true_graph(p,p);
	    NumericVector omega_true(p*(p-1)/2);
	    t = 0;

	    if(true_graph.isNull()){
	    	Rcout << "True graph is not available." << endl;
	    }
	    else{
	    	Rcout << "True graph is available." << endl;
	       rcpp_true_graph = as<NumericMatrix>(true_graph);
		   for (size_t i = 1; i < p; i++){
			   for (size_t j = 0; j < i; j++){
				   omega_true[t]=rcpp_true_graph(j,i);
				   t++;
			   }
		   }
	   }

	   // Get threshold
	   double upper = 2*sqrt(log(p));
	   NumericVector z(2000);
	   NumericVector numerator1(2000);
	   NumericVector denominator(2000);
	   for(double k = 0; k < 2000; k++){
		   z[k] = ((k+1)/2000) * upper;
		   numerator1[k] = 2 * (1-Rf_pnorm5(z[k], 0.0, 1.0, 1, 0)) *p*(p-1)/2;
		   denominator[k] = sum(abs(upper_stat)>=z[k]);
		   if(denominator[k] < 1){
			   denominator[k] = 1;
		   }
	   }
	   NumericVector FDP = numerator1/denominator;

	   NumericVector position(alpha_size);
	   for(size_t i = 0; i < alpha_size; i++){
		   for (double k = 0; k < 2000; k++){
			   if(FDP[k] <= rcpp_level[i]){		// threshold for FDP <= pre-specified level
			   			threshold[i] = z[k];
			   			position[i] = k;
			   			break;
			   }
			   else{
				   if(k == 1999){
					   threshold[i] = z[k];
					   position[i] = k;
				   }
			   }
		   }
	   }


	   if(true_graph.isNull()){
		   	   for(size_t i = 0; i < alpha_size; i++){
		   		   FDR[i] = FDP[position[i]];
		   	   }
		   	   for(size_t k = 0; k < alpha_size; k++){
		   		   NumericMatrix decision(p,p);
		   		   for(size_t i = 0; i < p-1; i++){
		   			   for(size_t j = i+1; j < p; j++){
		   				   if(abs(z_score(i,j)) >= threshold[k]){
		   					   decision(i,j) = 1;
		   				   }
		   				   decision(j,i) = decision(i,j);
		   			   }
		   		   }
		   		   global_decision[k] = decision;
		   	   }
	   	   }
	   	   else{
	   		   for(size_t i = 0; i < alpha_size; i++){
	   			   double FDR_denominator = sum(abs(upper_stat)>= threshold[i]);
	   			   if(FDR_denominator < 1){
	   				   	   FDR_denominator = 1;
	   			   }

	   			   FDR[i] = sum((abs(upper_stat) >= threshold[i]) & (omega_true == 0)) / FDR_denominator;

	   			   double Power_denominator = sum(omega_true!=0);
	   			   if(Power_denominator < 1){
	   				   Power_denominator = 1;
	   			   }
	   			   Power[i] = sum((abs(upper_stat) >= threshold[i]) & (omega_true != 0)) / Power_denominator;
	   		   }
		   	   for(size_t k = 0; k < alpha_size; k++){
		   		   NumericMatrix decision(p,p);
		   		   for(size_t i = 0; i < p-1; i++){
		   			   for(size_t j = i+1; j < p; j++){
		   				   if(abs(z_score(i,j)) >= threshold[k]){
		   					   decision(i,j) = 1;
		   				   }
		   				   decision(j,i) = decision(i,j);
		   			   }
		   		   }
		   		   global_decision[k] = decision;
		   	   }
	   	   }
  }

  if(global == false){
	  return List::create(_["precision"] = Omega, _["partialCor"] = partialCor, _["CI_low_precision"] = CI_low_Omega, _["CI_high_precision"] = CI_high_Omega, _["z_score_precision"] = z_score, _["p_precision"] = p_Omega);
  }
  else{
	  if(true_graph.isNull()){
		  return List::create(_["precision"] = Omega, _["partialCor"] = partialCor, _["CI_low_precision"] = CI_low_Omega, _["CI_high_precision"] = CI_high_Omega, _["z_score_precision"] = z_score, _["p_precision"] = p_Omega, _["threshold"] = threshold, _["FDR"] = FDR, _["global_decision"] = global_decision);
	  }
	  else{
		  return List::create(_["precision"] = Omega, _["partialCor"] = partialCor, _["CI_low_precision"] = CI_low_Omega, _["CI_high_precision"] = CI_high_Omega, _["z_score_precision"] = z_score, _["p_precision"] = p_Omega, _["threshold"] = threshold, _["FDR"] = FDR, _["power"] = Power, _["global_decision"] = global_decision);
	  }

  }
}

  if(rcpp_method[0] == "D-S_GL"){
	  NumericVector FDR(alpha_size);
	  NumericVector Power(alpha_size);
	  NumericVector threshold(alpha_size);
	  List global_decision(alpha_size);

  //penalty parameter for L1 norm of betas
  double rcpp_lambda;
  if(lambda.isNull()){
    rcpp_lambda = sqrt(log(p)/n);
    Rcout << "Use default lambda = sqrt(log(p)/n)" << endl;
    Rcout << "In this case, lambda = " << rcpp_lambda << endl;
  }
  else{
	rcpp_lambda = as<double>(lambda);
	Rcout << "In this case, lambda = " << rcpp_lambda << endl;
  }

  // Center matrix
  Rcout << "Center each column." << endl;
  NumericMatrix x_cent(n,p);
  for(size_t i = 0; i < p; i++){
    x_cent(_,i) = x(_,i) - mean(x(_,i));
  }


  // Empirical Sample Covariance Matrix
  Rcout << "Pre-calculate inner product matrixes." << endl;
  NumericMatrix S(p,p);
  for(size_t i = 0; i < p; i++){
	  for(size_t j = 0; j < p; j++){
		  S(i,j) = inner_product(x_cent(_,i).begin(), x_cent(_,i).end(), x_cent(_,j).begin(), 0.0)/n;
	  }
  }


  Rcout << "Calculate graphical Lasso." << endl;
  Environment pkg = Environment::namespace_env("glasso");
  Function f = pkg["glasso"];

  List fit = f(S, Named("rho") = rcpp_lambda, Named("penalize.diagonal") = false);

  NumericMatrix Omega = fit["wi"];

  // Symmetrize the Omega
  for(size_t i = 0; i < p-1; i++){
	  for(size_t j = i+1; j < p; j++){
		  Omega(i,j) = (Omega(i,j) + Omega(j,i))/2;
		  Omega(j,i) = Omega(i,j);
	  }
  }

  // Bias correction
  NumericMatrix inner1(p,p);
  NumericMatrix inner2(p,p);
  NumericMatrix Omega_correct(p,p);
  NumericMatrix partialCor(p,p);
  NumericMatrix CI_low_Omega_correct(p,p);
  NumericMatrix CI_high_Omega_correct(p,p);
  NumericMatrix z_score(p,p);
  NumericMatrix p_Omega_correct(p,p);

  for(size_t j = 0; j < p; j++){
	  size_t non_zero_count = 0;
	  NumericVector beta = Omega(_,j);
	  for (size_t k = 0; k < p; k++){
			 if(beta[k] != 0){
			     non_zero_count++;
			 }
	  }
	  NumericMatrix S_nonzero(p,non_zero_count);
	  NumericVector beta_nonzero(non_zero_count);

	  size_t t = 0;
	  for (size_t k = 0; k < p; k++){
		    if(beta[k] != 0){
		    	S_nonzero(_,t) = S(_,k);
		    	beta_nonzero[t] = beta[k];
		    	t++;
		    }
	  }

	  for(size_t i = 0; i < p; i++){
	     inner1(i,j) = inner_product(S_nonzero(i,_).begin(), S_nonzero(i,_).end(), beta_nonzero.begin(), 0.0);
	  }
  }

  for(size_t i = 0; i < p; i++){
	  size_t non_zero_count = 0;
	  NumericVector beta = Omega(i,_);
	  for (size_t k = 0; k < p; k++){
		  if(beta[k] != 0){
			  non_zero_count++;
		  }
	  }
	  NumericMatrix inner1_nonzero(non_zero_count, p);
	  NumericVector beta_nonzero(non_zero_count);

	  size_t t = 0;
		  for (size_t k = 0; k < p; k++){
			    if(beta[k] != 0){
			    	inner1_nonzero(t,_) = inner1(k,_);
			    	beta_nonzero[t] = beta[k];
			    	t++;
			    }
		  }


	  for(size_t j = 0; j < p; j++){
		  inner2(i,j) = inner_product(beta_nonzero.begin(), beta_nonzero.end(), inner1_nonzero(_,j).begin(),0.0);
		  Omega_correct(i,j) = 2*Omega(i,j) - inner2(i,j);
	  }
  }

  for(size_t i = 0; i < p-1; i++){
	  for(size_t j =i+1; j < p; j++){
		  partialCor(i,j) = -Omega_correct(i,j) / sqrt(Omega_correct(i,i)*Omega_correct(j,j));
		  double std_new = sqrt((Omega(i,i)*Omega(j,j)+pow(Omega(i,j),2))/n);
		  z_score(i,j) = Omega_correct(i,j)/std_new;
		  double z_95CI = 1.96;
		  CI_low_Omega_correct(i,j) = Omega_correct(i,j) - z_95CI * std_new;
		  CI_high_Omega_correct(i,j) = Omega_correct(i,j) + z_95CI * std_new;
		  double abs_zscore = abs(z_score(i,j));
		  p_Omega_correct(i,j) = 2 * Rf_pnorm5(-abs_zscore, 0.0, 1.0, 1, 0);

		  partialCor(j,i) = partialCor(i,j);
		  z_score(j,i) = z_score(i,j);
		  CI_low_Omega_correct(j,i) = CI_low_Omega_correct(i,j);
		  CI_high_Omega_correct(j,i) = CI_high_Omega_correct(i,j);
		  p_Omega_correct(j,i) = p_Omega_correct(i,j);
	  }
  }


  if(global == true){

  	    Rcout <<"Perform global inference." << endl;
	  	Rcout << "Use pre-specified level(s): " << rcpp_level << endl;

  	    //Take out the upper triangular of test statistic
  	    size_t t = 0;
  	    NumericVector upper_stat(p*(p-1)/2);
  	    for (size_t i = 1; i < p; i++){
  	    	for (size_t j = 0; j < i; j++){
  	    		upper_stat[t]=z_score(j,i);
  	    		t++;
  	    	}
  	    }

  	   //Take out the upper triangular of true graph if available
  	    NumericMatrix rcpp_true_graph(p,p);
  	    NumericVector omega_true(p*(p-1)/2);
  	    t = 0;

  	    if(true_graph.isNull()){
  	    	Rcout << "True graph is not available." << endl;
  	    }
  	    else{
  	    	Rcout << "True graph is available." << endl;
  	       rcpp_true_graph = as<NumericMatrix>(true_graph);
  		   for (size_t i = 1; i < p; i++){
  			   for (size_t j = 0; j < i; j++){
  				   omega_true[t]=rcpp_true_graph(j,i);
  				   t++;
  			   }
  		   }
  	   }

  	   // Get threshold
  	   double upper = 2*sqrt(log(p));
  	   NumericVector z(2000);
  	   NumericVector numerator1(2000);
  	   NumericVector denominator(2000);
  	   for(double k = 0; k < 2000; k++){
  		   z[k] = ((k+1)/2000) * upper;
  		   numerator1[k] = 2 * (1-Rf_pnorm5(z[k], 0.0, 1.0, 1, 0)) *p*(p-1)/2;
  		   denominator[k] = sum(abs(upper_stat)>=z[k]);
  		   if(denominator[k] < 1){
  			   denominator[k] = 1;
  		   }
  	   }
  	   NumericVector FDP = numerator1/denominator;

  	   NumericVector position(alpha_size);
	   for(size_t i = 0; i < alpha_size; i++){
		   for (double k = 0; k < 2000; k++){
			   if(FDP[k] <= rcpp_level[i]){		// threshold for FDP <= pre-specified level
			   			threshold[i] = z[k];
			   			position[i] = k;
			   			break;
			   }
			   else{
				   if(k == 1999){
					   threshold[i] = z[k];
					   position[i] = k;
				   }
			   }
		   }
	   }

  	 if(true_graph.isNull()){
	   	   	   for(size_t i = 0; i < alpha_size; i++){
	   	   		   FDR[i] = FDP[position[i]];
	   	   	   }
		   	   for(size_t k = 0; k < alpha_size; k++){
		   		   NumericMatrix decision(p,p);
		   		   for(size_t i = 0; i < p-1; i++){
		   			   for(size_t j = i+1; j < p; j++){
		   				   if(abs(z_score(i,j)) >= threshold[k]){
		   					   decision(i,j) = 1;
		   				   }
		   				   decision(j,i) = decision(i,j);
		   			   }
		   		   }
		   		   global_decision[k] = decision;
		   	   }
  	 }
  	 	   	   else{
  		   		   for(size_t i = 0; i < alpha_size; i++){
  		   			   double FDR_denominator = sum(abs(upper_stat)>= threshold[i]);
  		   			   if(FDR_denominator < 1){
  		   				   	   FDR_denominator = 1;
  		   			   }

  		   			   FDR[i] = sum((abs(upper_stat) >= threshold[i]) & (omega_true == 0)) / FDR_denominator;

  		   			   double Power_denominator = sum(omega_true!=0);
  		   			   if(Power_denominator < 1){
  		   				   Power_denominator = 1;
  		   			   }
  		   			   Power[i] = sum((abs(upper_stat) >= threshold[i]) & (omega_true != 0)) / Power_denominator;
  		   		   }
  			   	   for(size_t k = 0; k < alpha_size; k++){
  			   		   NumericMatrix decision(p,p);
  			   		   for(size_t i = 0; i < p-1; i++){
  			   			   for(size_t j = i+1; j < p; j++){
  			   				   if(abs(z_score(i,j)) >= threshold[k]){
  			   					   decision(i,j) = 1;
  			   				   }
  			   				   decision(j,i) = decision(i,j);
  			   			   }
  			   		   }
  			   		   global_decision[k] = decision;
  			   	   }
  	 	   	   }
  }

  if(global == false){
	  return List::create(_["precision"] = Omega_correct, _["partialCor"] = partialCor, _["CI_low_precision"] = CI_low_Omega_correct, _["CI_high_precision"] = CI_high_Omega_correct, _["z_score_precision"] = z_score, _["p_precision"] = p_Omega_correct);
  }
  else{
	  if(true_graph.isNull()){
		  return List::create(_["precision"] = Omega_correct, _["partialCor"] = partialCor, _["CI_low_precision"] = CI_low_Omega_correct, _["CI_high_precision"] = CI_high_Omega_correct, _["z_score_precision"] = z_score, _["p_precision"] = p_Omega_correct, _["threshold"] = threshold, _["FDR"] = FDR, _["global_decision"] = global_decision);
	  }
	  else{
		  return List::create(_["precision"] = Omega_correct, _["partialCor"] = partialCor, _["CI_low_precision"] = CI_low_Omega_correct, _["CI_high_precision"] = CI_high_Omega_correct, _["z_score_precision"] = z_score, _["p_precision"] = p_Omega_correct, _["threshold"] = threshold, _["FDR"] = FDR, _["power"] = Power, _["global_decision"] = global_decision);
	  }
  }
 }


  if(rcpp_method[0] == "GFC_SL"){
	  NumericMatrix T_value(p,p);
	  NumericMatrix test_statistic(p,p);
	  NumericVector FDR(alpha_size);
	  NumericVector Power(alpha_size);
	  NumericVector threshold(alpha_size);
	  List global_decision(alpha_size);

	  //penalty parameter for L1 norm of betas
	  double rcpp_lambda;
	  if(lambda.isNull()){
	    rcpp_lambda = sqrt(2 * log(p/sqrt(n))/n);
	    Rcout << "Use default lambda = sqrt(2*log(p/sqrt(n))/n)" << endl;
	    Rcout << "In this case, lambda = " << rcpp_lambda << endl;
	  }
	  else{
		rcpp_lambda = as<double>(lambda);
		Rcout << "In this case, lambda = " << rcpp_lambda << endl;
	  }

	  // Center matrix
	  Rcout << "Center each column." << endl;
	  NumericMatrix x_cent(n,p);
	  for(size_t i = 0; i < p; i++){
	    x_cent(_,i) = x(_,i) - mean(x(_,i));
	  }

	  // Standardize matrix
	  Rcout << "Standardize each column." << endl;
	  NumericVector x_l2norm(p);
	  for(size_t i = 0; i < p; i++){
	    x_l2norm[i] = sqrt(sum(pow(x_cent(_,i),2)));
	    if(x_l2norm[i] == 0){
	      stop("some variables have 0 variance, please remove them at first.");
	   }
	  }

	  NumericMatrix x_stan(n,p);
	  for(size_t i = 0; i < p; i++){
	    x_stan(_,i) = x_cent(_,i) / x_l2norm[i] * sqrt(n);
	  }


	  // Pre-calculate inner product matrixes
	  Rcout << "Pre-calculate inner product matrixes." << endl;
	  NumericMatrix IP_YX(p,p);
	  for(size_t i = 0; i < p; i++){
	    for(size_t j = 0; j < p; j++){
	      IP_YX(i,j) = inner_product(x_cent(_,i).begin(), x_cent(_,i).end(), x_stan(_,j).begin(), 0.0);
	    }
	  }
	  NumericMatrix IP_XX(p,p);
	  for(size_t i = 0; i < p; i++){
	    IP_XX(i,_) = IP_YX(i,_) / x_l2norm[i] * sqrt(n);
	  }

	  // Calculate scaled Lasso for each variable
	  Rcout << "Calculate scaled Lasso for each variable." << endl;
	  //NumericMatrix beta_pre(p,p-1);
	  NumericMatrix beta_initial(p,p);
	  NumericMatrix epsi_pre(n,p);

	  for(size_t i = 0; i < p; i++){
	    if(i % 100 == 0){
	      Rcout <<"scaled Lasso for variable " << (i + 1) << endl;
	    }
	    double sigma = 1.0;
	    NumericVector beta(p-1);
	    NumericVector epsi(n);
	    double reg, diffBeta, sigma_tmp, diffSigma;

	    // Extract IP_XX[-i,-i] and IP_YX[i,-i]
	    NumericMatrix tmpxx(p-1,p-1);
	    NumericVector tmpyx(p-1);
	    size_t t = 0;
	    for (size_t k = 0; k < p; k++){
	      if (k != i){
	        tmpyx[t] = IP_YX(i,k);
	        size_t t1 = 0;
	        for (size_t l = 0; l < p; l++){ // l is index of row, k is index of col
	          if (l != i){
	            tmpxx(t1,t) = IP_XX(l,k);
	            t1++;
	          }
	        }
	        t++;
	      }
	    }

	    // Extract x_stan[,-i]
	    NumericMatrix tmp_x(n, p-1);
	    t=0;
	    for (size_t k = 0; k < p; k++){
	      if (k != i){
	        tmp_x(_,t) = x_stan(_,k);
	        t++;
	      }
	    }

	    // Extract x_l2norm[,-i]
	    NumericVector x_l2norm_new(p-1);
	    t=0;
	    for (size_t k = 0; k < p; k++){
	      if (k != i){
	    	x_l2norm_new[t] = x_l2norm[k];
	    	t++;
	      }
	    }

	    //Scaled_lasso
	    size_t iter = 0;
	    do {
	      //Update beta when fix sigma
	      reg = sigma * rcpp_lambda;
	      NumericVector beta_tmp = FastLasso_Rcpp(tmpyx, tmpxx, reg, n);
	      diffBeta = max(abs(beta_tmp - beta));
	      beta = beta_tmp;
	      //Update sigma when fix beta
	      NumericVector tmp_y(n);

	      size_t non_zero_count = 0;
	            for (size_t k = 0; k < p-1; k++){
	          	  if(beta[k] != 0){
	          		  non_zero_count++;
	          	  }
	            }

	            NumericMatrix tmp_x_nonzero(n,non_zero_count);
	            NumericVector beta_nonzero(non_zero_count);

	            size_t t = 0;
	            for (size_t k = 0; k < p-1; k++){
	          	  if(beta[k] != 0){
	          		  tmp_x_nonzero(_,t) = tmp_x(_,k);
	          		  beta_nonzero[t] = beta[k];
	          		  t++;
	          	  }
	            }

	            for (size_t k = 0; k < n; k++){
	              tmp_y[k] = inner_product(tmp_x_nonzero(k,_).begin(), tmp_x_nonzero(k,_).end(), beta_nonzero.begin(), 0.0);
	            } //Multiplication of x.stan[,-i] and beta

	      epsi = x_cent(_,i) - tmp_y;
	      sigma_tmp = sqrt(sum(pow(epsi,2))/n);
	      diffSigma = abs(sigma_tmp - sigma);
	      sigma = sigma_tmp;
	      iter++;
	    }
	    while((diffSigma >= tol || diffBeta >= tol) && iter < 10); //Stop iteration

	    epsi_pre(_,i) = epsi;

	    for (size_t k = 0; k < p; k++){
	    	if(i < k){
	    		beta_initial(i,k) = beta[k-1]/x_l2norm_new[k-1]*sqrt(n);
	    	}
	    	if(i > k){
	    		beta_initial(i,k) = beta[k]/x_l2norm_new[k]*sqrt(n);
	    	}
	    }

	  }

	  NumericMatrix inner_epsi(p,p);
	  for (size_t i = 0; i < p; i++){
		  for (size_t j = 0; j < p; j++){
			  inner_epsi(i,j) = inner_product(epsi_pre(_,i).begin(), epsi_pre(_,i).end(), epsi_pre(_,j).begin(), 0.0);
		  }
	  }

	    //Bias Correction

	    for (size_t i = 0; i < p-1; i++){
	    	for(size_t j=i+1; j < p; j++){

	            T_value(i,j)=(inner_epsi(i,j)+inner_epsi(i,i)*beta_initial(j,i)+inner_epsi(j,j)*beta_initial(i,j))/n;
	            double std_new = sqrt((inner_epsi(i,i)/n)*(inner_epsi(j,j)/n)/n);
	            test_statistic(i,j)=T_value(i,j)/std_new;

	            T_value(j,i)=T_value(i,j);
	            test_statistic(j,i)=test_statistic(i,j);

	   	}
	  }

	    Rcout <<"Perform global inference." << endl;
	    Rcout << "Use pre-specified level(s): " << rcpp_level << endl;

	   //Take out the upper triangular of test statistic
	    size_t t = 0;
	    NumericVector upper_stat(p*(p-1)/2);
	    for (size_t i = 1; i < p; i++){
	    	for (size_t j = 0; j < i; j++){
	    		upper_stat[t]=test_statistic(j,i);
	    		t++;
	    	}
	    }

	   //Take out the upper triangular of true graph if available
	   	NumericMatrix rcpp_true_graph(p,p);
	   	NumericVector omega_true(p*(p-1)/2);
	   	t = 0;
	    if(true_graph.isNull()){
	    	Rcout << "True graph is not available." << endl;
	    }
	    else{
	    	Rcout << "True graph is available." << endl;
	       rcpp_true_graph = as<NumericMatrix>(true_graph);
		   for (size_t i = 1; i < p; i++){
			   for (size_t j = 0; j < i; j++){
				   omega_true[t]=rcpp_true_graph(j,i);
				   t++;
			   }
		   }
	   }


	   // Get threshold
	   double upper = 2*sqrt(log(p));
	   NumericVector z(2000);
	   NumericVector numerator1(2000);
	   NumericVector denominator(2000);
	   for(double k = 0; k < 2000; k++){
		   z[k] = ((k+1)/2000) * upper;
		   numerator1[k] = 2 * (1-Rf_pnorm5(z[k], 0.0, 1.0, 1, 0)) *p*(p-1)/2;
		   denominator[k] = sum(abs(upper_stat)>=z[k]);
		   if(denominator[k] < 1){
			   denominator[k] = 1;
		   }
	   }
	   NumericVector FDP = numerator1/denominator;

	   NumericVector position(alpha_size);
	   for(size_t i = 0; i < alpha_size; i++){
		   for (double k = 0; k < 2000; k++){
			   if(FDP[k] <= rcpp_level[i]){		// threshold for FDP <= pre-specified level
			   			threshold[i] = z[k];
			   			position[i] = k;
			   			break;
			   }
			   else{
				   if(k == 1999){
					   threshold[i] = z[k];
					   position[i] = k;
				   }
			   }
		   }
	   }


	   if(true_graph.isNull()){
   	   	   for(size_t i = 0; i < alpha_size; i++){
	   	   		   FDR[i] = FDP[position[i]];
	   	   	   }
	   	   for(size_t k = 0; k < alpha_size; k++){
	   		   NumericMatrix decision(p,p);
	   		   for(size_t i = 0; i < p-1; i++){
	   			   for(size_t j = i+1; j < p; j++){
	   				   if(abs(test_statistic(i,j)) >= threshold[k]){
	   					   decision(i,j) = 1;
	   				   }
	   				   decision(j,i) = decision(i,j);
	   			   }
	   		   }
	   		   global_decision[k] = decision;
	   	   }
	   }
	   else{

	   		   for(size_t i = 0; i < alpha_size; i++){
	   			   double FDR_denominator = sum(abs(upper_stat)>= threshold[i]);
	   			   if(FDR_denominator < 1){
	   				   	   FDR_denominator = 1;
	   			   }

	   			   FDR[i] = sum((abs(upper_stat) >= threshold[i]) & (omega_true == 0)) / FDR_denominator;

	   			   double Power_denominator = sum(omega_true!=0);
	   			   if(Power_denominator < 1){
	   				   Power_denominator = 1;
	   			   }
	   			   Power[i] = sum((abs(upper_stat) >= threshold[i]) & (omega_true != 0)) / Power_denominator;
	   		   }
		   	   for(size_t k = 0; k < alpha_size; k++){
		   		   NumericMatrix decision(p,p);
		   		   for(size_t i = 0; i < p-1; i++){
		   			   for(size_t j = i+1; j < p; j++){
		   				   if(abs(test_statistic(i,j)) >= threshold[k]){
		   					   decision(i,j) = 1;
		   				   }
		   				   decision(j,i) = decision(i,j);
		   			   }
		   		   }
		   		   global_decision[k] = decision;
		   	   }
	   }

	   if(true_graph.isNull()){
		   return List::create(_["T_stat"] = test_statistic, _["threshold"] = threshold, _["FDR"] = FDR, _["global_decision"] = global_decision);
	   }
	   else{
		   return List::create(_["T_stat"] = test_statistic, _["threshold"] = threshold, _["FDR"] = FDR, _["power"] = Power, _["global_decision"] = global_decision);
	   }
  }

  if(rcpp_method[0] == "GFC_L"){
	  NumericVector FDR(alpha_size);
	  NumericVector Power(alpha_size);
	  NumericVector threshold(alpha_size);
	  List global_decision(alpha_size);

	  // Center matrix
	  Rcout << "Center each column." << endl;
	  NumericMatrix x_cent(n,p);
	  for(size_t i = 0; i < p; i++){
	    x_cent(_,i) = x(_,i) - mean(x(_,i));
	  }

	  // Standardize matrix
	  Rcout << "Standardize each column." << endl;
	  NumericVector x_l2norm(p);
	  for(size_t i = 0; i < p; i++){
	    x_l2norm[i] = sqrt(sum(pow(x_cent(_,i),2)));
	    if(x_l2norm[i] == 0){
	      stop("some variables have 0 variance, please remove them at first.");
	   }
	  }

	  NumericMatrix x_stan(n,p);
	  for(size_t i = 0; i < p; i++){
	    x_stan(_,i) = x_cent(_,i) / x_l2norm[i] * sqrt(n);
	  }



	  // Pre-calculate inner product matrixes
	  Rcout << "Pre-calculate inner product matrixes." << endl;
	  NumericMatrix IP_YX(p,p);
	  for(size_t i = 0; i < p; i++){
	    for(size_t j = 0; j < p; j++){
	      IP_YX(i,j) = inner_product(x_cent(_,i).begin(), x_cent(_,i).end(), x_stan(_,j).begin(), 0.0);
	    }
	  }
	  NumericMatrix IP_XX(p,p);
	  NumericMatrix IP_YX_origin(p,p);
	  for(size_t i = 0; i < p; i++){
	    IP_XX(i,_) = IP_YX(i,_) / x_l2norm[i] * sqrt(n);
	    IP_YX_origin(_,i) = IP_YX(_,i) * x_l2norm[i] / sqrt(n);
	  }


	  // Calculate Lasso for each variable
	  NumericMatrix beta_initial(p,p);
	  NumericMatrix epsi_pre(n,p);

	  double rcpp_ndelta;
	  if(ndelta.isNull()){
		  rcpp_ndelta = 40;
		  Rcout << "Use default number of delta = 40" << endl;
	  }
	  else{
		  rcpp_ndelta = as<double>(ndelta);
		  Rcout << "In this case, number of delta = " << rcpp_ndelta << endl;
	  }


	  // Define a sequence of delta from 0 to 2
	  NumericVector delta(rcpp_ndelta);
	  for(double m = 0; m < rcpp_ndelta; m++){
		  delta[m] = (m+1)/(rcpp_ndelta/2);
	  }

	  delta = sortDescending(delta);


	  NumericVector error_cum(rcpp_ndelta);
	  List T(rcpp_ndelta);
	  List Statistic(rcpp_ndelta);
	  List upper_stat_tmp(rcpp_ndelta);

	    Rcout <<"Perform global inference." << endl;
	    Rcout << "Use pre-specified level(s): " << rcpp_level << endl;


	   //Take out the upper triangular of true graph if available
	   	NumericMatrix rcpp_true_graph(p,p);
	   	NumericVector omega_true(p*(p-1)/2);
	   	size_t t = 0;
	    if(true_graph.isNull()){
	    	Rcout << "True graph is not available." << endl;
	    }
	    else{
	    	Rcout << "True graph is available." << endl;
	       rcpp_true_graph = as<NumericMatrix>(true_graph);
		   for (size_t i = 1; i < p; i++){
			   for (size_t j = 0; j < i; j++){
				   omega_true[t]=rcpp_true_graph(j,i);
				   t++;
			   }
		   }
	   }


	  NumericMatrix tmp_beta(p-1,p);
	  NumericMatrix tmp_gc(p-1,p);
	  size_t count = 0;

	  Rcout << "Calculate Lasso of each variable with tuning parameters under each delta." << endl;
	  Rcout << "Record test statistics under each delta." << endl;
	  for(size_t m = 0; m < rcpp_ndelta; m++){
		  NumericMatrix T_value(p,p);
		  NumericMatrix test_statistic(p,p);

		  count++;

		  Rcout <<"delta " << (count) << endl;

		  for(size_t i = 0; i < p; i++){

			  double lambda = delta[m]*sqrt(IP_YX_origin(i,i)/n*(log(p)/n));

			  if(i % 100 == 0){
				  Rcout <<"Lasso for variable " << (i + 1) << endl;
			  }

			  NumericVector epsi(n);

			  // Extract IP_XX[-i,-i] and IP_YX[i,-i]
			  NumericMatrix tmpxx(p-1,p-1);
			  NumericVector tmpyx(p-1);
			  size_t t = 0;
			  for (size_t k = 0; k < p; k++){
				  if (k != i){
					  tmpyx[t] = IP_YX(i,k);
					  size_t t1 = 0;
					  for (size_t l = 0; l < p; l++){ // l is index of row, k is index of col
						  if (l != i){
							  tmpxx(t1,t) = IP_XX(l,k);
							  t1++;
						  }
					  }
					  t++;
				  }
			  }

			  // Extract x_stan[,-i]
			  NumericMatrix tmp_x(n, p-1);
			  t=0;
			  for (size_t k = 0; k < p; k++){
				  if (k != i){
					  tmp_x(_,t) = x_stan(_,k);
					  t++;
				  }
			  }


			  // Extract x_l2norm[,-i]
			  NumericVector x_l2norm_new(p-1);
			  t=0;
			  for (size_t k = 0; k < p; k++){
				  if (k != i){
					  x_l2norm_new[t] = x_l2norm[k];
					  t++;
				  }
			  }

			  // Lasso
			  NumericVector beta_column = tmp_beta(_,i);
			  NumericVector gc_column = tmp_gc(_,i);

			  List beta_result = FastLasso_Rcpp_v2(beta_column, gc_column, tmpyx, tmpxx, lambda, n);

			  NumericVector beta = beta_result["beta"];

			  NumericVector gc = beta_result["gc"];

			  tmp_beta(_,i) = beta;
			  tmp_gc(_,i) = gc;

			  //Multiplication of x.stan[,-i] and beta
			  NumericVector tmp_y(n);
			  size_t non_zero_count = 0;
			  for (size_t k = 0; k < p-1; k++){
			  	          	  if(beta[k] != 0){
			  	          		  non_zero_count++;
			  	          	  }
			  	            }

			  	            NumericMatrix tmp_x_nonzero(n,non_zero_count);
			  	            NumericVector beta_nonzero(non_zero_count);

			  	            t = 0;
			  	            for (size_t k = 0; k < p-1; k++){
			  	          	  if(beta[k] != 0){
			  	          		  tmp_x_nonzero(_,t) = tmp_x(_,k);
			  	          		  beta_nonzero[t] = beta[k];
			  	          		  t++;
			  	          	  }
			  	            }

			  for (size_t k = 0; k < n; k++ ){
				  tmp_y[k] = inner_product(tmp_x_nonzero(k,_).begin(), tmp_x_nonzero(k,_).end(), beta_nonzero.begin(), 0.0);
			  }

			  //Record residual
			  epsi = x_cent(_,i) - tmp_y;

			  epsi_pre(_,i) = epsi;

			  for (size_t k = 0; k < p; k++){
				  if(i < k){
					  beta_initial(i,k) = beta[k-1]/x_l2norm_new[k-1]*sqrt(n);
				  }
				  if(i > k){
					  beta_initial(i,k) = beta[k]/x_l2norm_new[k]*sqrt(n);
				  }
			  }
		  }


		  NumericMatrix inner_epsi(p,p);
		  for (size_t i = 0; i < p; i++){
			  for (size_t j = 0; j < p; j++){
				  inner_epsi(i,j) = inner_product(epsi_pre(_,i).begin(), epsi_pre(_,i).end(), epsi_pre(_,j).begin(), 0.0);
			  }
		  }

		  //Bias Correction

		  for (size_t i = 0; i < p-1; i++){
			  for(size_t j=i+1; j < p; j++){

	            T_value(i,j)=(inner_epsi(i,j)+inner_epsi(i,i)*beta_initial(j,i)+inner_epsi(j,j)*beta_initial(i,j))/n;
	            double std_new = sqrt((inner_epsi(i,i)/n)*(inner_epsi(j,j)/n)/n);
	            test_statistic(i,j)=T_value(i,j)/std_new;

	            T_value(j,i)=T_value(i,j);
	            test_statistic(j,i)=test_statistic(i,j);

			  }
		  }

		  T[m] = T_value;
		  Statistic[m] = test_statistic;

		  //Take out the upper triangular of test statistic
		  size_t t = 0;
		  NumericVector upper_stat(p*(p-1)/2);
		  for (size_t i = 1; i < p; i++){
			  for (size_t j = 0; j < i; j++){
				  upper_stat[t]=test_statistic(j,i);
				  t++;
			  }
		  }


		  //threshold for choosing delta
		  double error_1 = 0;
		  for (double k = 3; k < 10; k++){
			  double normal_tail = Rf_qnorm5(1-k/20, 0.0, 1.0, 1, 0);
			  double error = sum(abs(upper_stat)>=normal_tail)/(p*(p-1)*k/20)-1;
			  error_1 = error_1 + pow(error,2);
		  }

		  upper_stat_tmp[m] = upper_stat;
		  error_cum[m] = error_1;
	  }

	  	 //Reverse the error_cum
	  	 reverse(error_cum.begin(), error_cum.end());


	  	 //Return the minimum position for error

	  	 Rcout << "Choose delta for FDR control." << endl;
		double position_select = min_position(error_cum);

		NumericMatrix T_value = T[rcpp_ndelta-position_select-1];
		NumericMatrix test_stat = Statistic[rcpp_ndelta-position_select-1];
		NumericVector upper_stat_select = upper_stat_tmp[rcpp_ndelta-position_select-1];


		  // Get threshold
		  double upper = 2*sqrt(log(p));
		  NumericVector z(2000);
		  NumericVector numerator1(2000);
		  NumericVector denominator(2000);
		  for(double k = 0; k < 2000; k++){
			   z[k] = ((k+1)/2000) * upper;
			   numerator1[k] = 2 * (1-Rf_pnorm5(z[k], 0.0, 1.0, 1, 0)) *p*(p-1)/2;
			   denominator[k] = sum(abs(upper_stat_select)>=z[k]);
			   if(denominator[k] < 1){
				   denominator[k] = 1;
			   }
		  }
		  NumericVector FDP = numerator1/denominator;


		   NumericVector position(alpha_size);

		   for(size_t i = 0; i < alpha_size; i++){
			   for (double k = 0; k < 2000 ; k++){
				   if(FDP[k] <= rcpp_level[i]){      // threshold for FDP <= pre-specified level
					   	   threshold[i] = z[k];
					   	   position[i] = k;
					   	   break;
				   }
				   else{
					   if(k==1999){
						   threshold[i] = z[k];
					   	   position[i] = k;
					   }
				   }
			   }
		   }

		   if(true_graph.isNull()){
			   for(size_t i = 0; i < alpha_size; i++){
			   	   		   FDR[i] = FDP[position[i]];
			   }
		   	   for(size_t k = 0; k < alpha_size; k++){
		   		   NumericMatrix decision(p,p);
		   		   for(size_t i = 0; i < p-1; i++){
		   			   for(size_t j = i+1; j < p; j++){
		   				   if(abs(test_stat(i,j)) >= threshold[k]){
		   					   decision(i,j) = 1;
		   				   }
		   				   decision(j,i) = decision(i,j);
		   			   }
		   		   }
		   		   global_decision[k] = decision;
		   	   }
		   }
		   else{

	   		   for(size_t i = 0; i < alpha_size; i++){
	   			   double FDR_denominator = sum(abs(upper_stat_select)>= threshold[i]);
	   			   if(FDR_denominator < 1){
	   				   	   FDR_denominator = 1;
	   			   }

	   			   FDR[i] = sum((abs(upper_stat_select) >= threshold[i]) & (omega_true == 0)) / FDR_denominator;

	   			   double Power_denominator = sum(omega_true!=0);
	   			   if(Power_denominator < 1){
	   				   Power_denominator = 1;
	   			   }
	   			   Power[i] = sum((abs(upper_stat_select) >= threshold[i]) & (omega_true != 0)) / Power_denominator;
	   		   }
		   	   for(size_t k = 0; k < alpha_size; k++){
		   		   NumericMatrix decision(p,p);
		   		   for(size_t i = 0; i < p-1; i++){
		   			   for(size_t j = i+1; j < p; j++){
		   				   if(abs(test_stat(i,j)) >= threshold[k]){
		   					   decision(i,j) = 1;
		   				   }
		   				   decision(j,i) = decision(i,j);
		   			   }
		   		   }
		   		   global_decision[k] = decision;
		   	   }
		   }



	  if(true_graph.isNull()){
		  return List::create(_["T_stat"] = test_stat, _["threshold"] = threshold, _["FDR"] = FDR, _["global_decision"] = global_decision);
	  }
	  else{
		  return List::create(_["T_stat"] = test_stat, _["threshold"] = threshold, _["FDR"] = FDR, _["power"] = Power, _["global_decision"] = global_decision);
	   }
  }

  if(rcpp_method[0] == "B_NW_SL"){
	  //penalty parameter for L1 norm of betas
	  double rcpp_lambda;
	  if(lambda.isNull()){
		rcpp_lambda = sqrt(2 * log(p/sqrt(n))/n);
	    Rcout << "Use default lambda = sqrt(2*log(p/sqrt(n))/n)" << endl;
	    Rcout << "In this case, lambda = " << rcpp_lambda << endl;
	  }
	  else{
		 rcpp_lambda = as<double>(lambda);
		 Rcout << "In this case, lambda = " << rcpp_lambda << endl;
	  }

	  NumericMatrix precision(p,p);
	  NumericMatrix p_precision(p,p);
	  NumericMatrix partialCor(p,p);
	  NumericMatrix p_partialCor(p,p);
	  NumericMatrix CI_low_parCor(p,p);
	  NumericMatrix CI_high_parCor(p,p);
	  NumericMatrix CI_low_precision(p,p);
	  NumericMatrix CI_high_precision(p,p);
	  NumericMatrix z_score_precision(p,p);
	  NumericMatrix z_score_partialCor(p,p);
	  NumericVector FDR(alpha_size);
	  NumericVector Power(alpha_size);
	  NumericVector threshold(alpha_size);
	  List global_decision(alpha_size);


	  // Center matrix
	  Rcout << "Center each column." << endl;
	  NumericMatrix x_cent(n,p);
	  for(size_t i = 0; i < p; i++){
	    x_cent(_,i) = x(_,i) - mean(x(_,i));
	  }

	  // Standardize matrix
	  Rcout << "Standardize each column." << endl;
	  NumericVector x_l2norm(p);
	  for(size_t i = 0; i < p; i++){
	    x_l2norm[i] = sqrt(sum(pow(x_cent(_,i),2)));
	    if(x_l2norm[i] == 0){
	      stop("some variables have 0 variance, please remove them at first.");
	    }
	  }
	  NumericMatrix x_stan(n,p);
	  for(size_t i = 0; i < p; i++){
	    x_stan(_,i) = x_cent(_,i) / x_l2norm[i] * sqrt(n);
	  }

	  // Pre-calculate inner product matrixes
	  Rcout << "Pre-calculate inner product matrixes." << endl;
	  NumericMatrix IP_YX(p,p);
	  for(size_t i = 0; i < p; i++){
	    for(size_t j = 0; j < p; j++){
	      IP_YX(i,j) = inner_product(x_cent(_,i).begin(), x_cent(_,i).end(), x_stan(_,j).begin(), 0.0);
	    }
	  }
	  NumericMatrix IP_XX(p,p);
	  for(size_t i = 0; i < p; i++){
	    IP_XX(i,_) = IP_YX(i,_) / x_l2norm[i] * sqrt(n);
	  }

	  // Pre-calculate scaled Lasso for each variable
	  Rcout << "Pre-calculate scaled Lasso for each variable." << endl;
	  NumericMatrix beta_pre(p,p);
	  NumericMatrix epsi_pre(n,p);

	  for(size_t i = 0; i < p; i++){
	    if(i % 100 == 0){
	      Rcout <<"Pre-Lasso for variable " << (i + 1) << endl;
	    }
	    double sigma = 1.0;
	    NumericVector beta(p-1);
	    NumericVector epsi(n);
	    double reg, diffBeta, sigma_tmp, diffSigma;

	    // Extract IP_XX[-i,-i] and IP_YX[i,-i]
	    NumericMatrix tmpxx(p-1,p-1);
	    NumericVector tmpyx(p-1);
	    size_t t = 0;
	    for (size_t k = 0; k < p; k++){
	      if (k != i){
	        tmpyx[t] = IP_YX(i,k);
	        size_t t1 = 0;
	        for (size_t l = 0; l < p; l++){ // l is index of row, k is index of col
	          if (l != i){
	            tmpxx(t1,t) = IP_XX(l,k);
	            t1++;
	          }
	        }
	        t++;
	      }
	    }
	    // Extract x_stan[,-i]
	    NumericMatrix tmp_x(n, p-1);
	    t=0;
	    for (size_t k = 0; k < p; k++){
	      if (k != i){
	        tmp_x(_,t) = x_stan(_,k);
	        t++;
	      }
	    }

	    //Scaled_lasso
	    size_t iter = 0;
	    do {
	      //Update beta when fix sigma
	      reg = sigma * rcpp_lambda;
	      NumericVector beta_tmp = FastLasso_Rcpp(tmpyx, tmpxx, reg, n);
	      diffBeta = max(abs(beta_tmp - beta));
	      beta = beta_tmp;
	      //Update sigma when fix beta
	      NumericVector tmp_y(n);
	      size_t non_zero_count = 0;
	            for (size_t k = 0; k < p-1; k++){
	          	  if(beta[k] != 0){
	          		  non_zero_count++;
	          	  }
	            }

	            NumericMatrix tmp_x_nonzero(n,non_zero_count);
	            NumericVector beta_nonzero(non_zero_count);

	            size_t t = 0;
	            for (size_t k = 0; k < p-1; k++){
	          	  if(beta[k] != 0){
	          		  tmp_x_nonzero(_,t) = tmp_x(_,k);
	          		  beta_nonzero[t] = beta[k];
	          		  t++;
	          	  }
	            }

	            for (size_t k = 0; k < n; k++){
	              tmp_y[k] = inner_product(tmp_x_nonzero(k,_).begin(), tmp_x_nonzero(k,_).end(), beta_nonzero.begin(), 0.0);
	            } //Multiplication of x.stan[,-i] and beta

	      epsi = x_cent(_,i) - tmp_y;
	      sigma_tmp = sqrt(sum(pow(epsi,2))/n);
	      diffSigma = abs(sigma_tmp - sigma);
	      sigma = sigma_tmp;
	      iter++;
	    }
	    while((diffSigma >= tol || diffBeta >= tol) && iter < 10); //Stop iteration
	    epsi_pre(_,i) = epsi;
	    t = 0;
	    for (size_t k = 0; k < p; k++){
	      if (k != i){
	        beta_pre(k,i) = beta[t];
	        t++;
	      }
	    }
	  }

	  // Fast pairwise GGM based on the pre-calculated scaled Lasso
	  Rcout << "Perform pairwise GGM." << endl;
	  for(size_t i = 0; i < (p-1); i++){
	    if(i % 100 == 0){
	      Rcout <<"Pair-Lasso for variable " << (i + 1) << endl;
	    }
	    for(size_t j = (i+1); j < p; j++){
	      NumericVector epsi_i(n);
	      NumericVector epsi_j(n);

	      // Solve scaled Lasso i without variable j
	      if(beta_pre(j,i) == 0){
	        epsi_i = epsi_pre(_,i);
	      } else{
	        double sigma = 1.0;
	        NumericVector beta(p-2);
	        NumericVector epsi(n);
	        double reg, diffBeta, sigma_tmp, diffSigma;

	        // Extract IP_XX[-c(i,j),-c(i,j)] and IP_YX[i,-c(i,j)]
	        NumericMatrix tmpxx(p-2,p-2);
	        NumericVector tmpyx(p-2);
	        size_t t = 0;
	        for (size_t k = 0; k < p; k++){
	          if (k != i && k != j){
	            tmpyx[t] = IP_YX(i,k);
	            size_t t1 = 0;
	            for (size_t l = 0; l < p; l++){
	              if ( l != i && l != j){
	                tmpxx(t1,t) = IP_XX(l,k);
	                t1++;
	              }
	            }
	            t++;
	          }
	        }
	        // Extract x_stan[,-c(i,j)]
	        NumericMatrix tmp_x(n, p-2);
	        t=0;
	        for (size_t k = 0; k < p; k++){
	          if (k != i && k != j){
	            tmp_x(_,t) = x_stan(_,k);
	            t++;
	          }
	        }

	        //Scaled_lasso
	        size_t iter = 0;
	        do {
	          //Update beta when fix sigma
	          reg = sigma * rcpp_lambda;
	          NumericVector beta_tmp = FastLasso_Rcpp(tmpyx, tmpxx, reg, n);
	          diffBeta = max(abs(beta_tmp - beta));
	          beta = beta_tmp;
	          //Update sigma when fix beta
	          NumericVector tmp_y(n);
	          size_t non_zero_count = 0;
	                for (size_t k = 0; k < p-2; k++){
	              	  if(beta[k] != 0){
	              		  non_zero_count++;
	              	  }
	                }

	                NumericMatrix tmp_x_nonzero(n,non_zero_count);
	                NumericVector beta_nonzero(non_zero_count);

	                size_t t = 0;
	                for (size_t k = 0; k < p-2; k++){
	              	  if(beta[k] != 0){
	              		  tmp_x_nonzero(_,t) = tmp_x(_,k);
	              		  beta_nonzero[t] = beta[k];
	              		  t++;
	              	  }
	                }

	                for (size_t k = 0; k < n; k++){
	                  tmp_y[k] = inner_product(tmp_x_nonzero(k,_).begin(), tmp_x_nonzero(k,_).end(), beta_nonzero.begin(), 0.0);
	                } //Multiplication of x.stan[,-i] and beta

	                epsi = x_cent(_,i) - tmp_y;
	          sigma_tmp = sqrt(sum(pow(epsi,2))/n);
	          diffSigma = abs(sigma_tmp - sigma);
	          sigma = sigma_tmp;
	          iter++;
	        }
	        while((diffSigma >= tol || diffBeta >= tol) && iter < 10);
	        epsi_i = epsi;
	      }

	      // Solve scaled Lasso j without variable i
	      if(beta_pre(i,j) == 0){
	        epsi_j = epsi_pre(_,j);
	      } else{
	        double sigma = 1.0;
	        NumericVector beta(p-2);
	        NumericVector epsi(n);
	        double reg, diffBeta, sigma_tmp, diffSigma;

	        NumericMatrix tmpxx(p-2,p-2);
	        NumericVector tmpyx(p-2);
	        size_t t = 0;
	        for (size_t k = 0; k < p; k++){
	          if ( k != i && k != j){
	            tmpyx[t] = IP_YX(j,k);
	            size_t t1 = 0;
	            for (size_t l = 0; l < p; l++){
	              if ( l != i && l != j){
	                tmpxx(t1,t) = IP_XX(l,k);
	                t1++;
	              }
	            }
	            t++;
	          }
	        }
	        NumericMatrix tmp_x(n, p-2);
	        t = 0;
	        for (size_t k = 0; k < p; k++){
	          if (k != i && k != j){
	            tmp_x(_,t) = x_stan(_,k);
	            t++;
	          }
	        }

	        //Scaled_lasso
	        size_t iter = 0;
	        do {
	          //Update beta when fix sigma
	          reg = sigma * rcpp_lambda;
	          NumericVector beta_tmp = FastLasso_Rcpp(tmpyx, tmpxx, reg, n);
	          diffBeta = max(abs(beta_tmp - beta));
	          beta = beta_tmp;
	          //Update sigma when fix beta
	          NumericVector tmp_y(n);
	          size_t non_zero_count = 0;
	                for (size_t k = 0; k < p-2; k++){
	              	  if(beta[k] != 0){
	              		  non_zero_count++;
	              	  }
	                }

	                NumericMatrix tmp_x_nonzero(n,non_zero_count);
	                NumericVector beta_nonzero(non_zero_count);

	                size_t t = 0;
	                for (size_t k = 0; k < p-2; k++){
	              	  if(beta[k] != 0){
	              		  tmp_x_nonzero(_,t) = tmp_x(_,k);
	              		  beta_nonzero[t] = beta[k];
	              		  t++;
	              	  }
	                }

	                for (size_t k = 0; k < n; k++){
	                  tmp_y[k] = inner_product(tmp_x_nonzero(k,_).begin(), tmp_x_nonzero(k,_).end(), beta_nonzero.begin(), 0.0);
	                } //Multiplication of x.stan[,-i] and beta
	          epsi = x_cent(_,j) - tmp_y;
	          sigma_tmp = sqrt(sum(pow(epsi,2))/n);
	          diffSigma = abs(sigma_tmp - sigma);
	          sigma = sigma_tmp;
	          iter++;
	        }
	        while((diffSigma >= tol || diffBeta >= tol) && iter < 10);
	        epsi_j = epsi;
	      }

	      // Precision, solve(t(epsi_ij)%*%epsi_ij/n)
	      NumericMatrix omega_tmp(2,2);
	      NumericMatrix omega(2,2);
	      omega_tmp(0,0) = inner_product(epsi_i.begin(), epsi_i.end(), epsi_i.begin(), 0.0)/n;
	      omega_tmp(1,1) = inner_product(epsi_j.begin(), epsi_j.end(), epsi_j.begin(), 0.0)/n;
	      omega_tmp(0,1) = inner_product(epsi_i.begin(), epsi_i.end(), epsi_j.begin(), 0.0)/n;
	      omega_tmp(1,0) = omega_tmp(0,1);
	      // Inverse matrix of omega_tmp
	      double tmp = omega_tmp(0,0) * omega_tmp(1,1) - omega_tmp(0,1) * omega_tmp(1,0);
	      omega(0,0) = omega_tmp(1,1)/tmp;
	      omega(0,1) = -omega_tmp(0,1)/tmp;
	      omega(1,0) = -omega_tmp(1,0)/tmp;
	      omega(1,1) = omega_tmp(0,0)/tmp;
	      precision(i,i) = omega(0,0);
	      precision(j,j) = omega(1,1);
	      precision(i,j) = omega(0,1);
	      precision(j,i) = omega(1,0);

	      // P-value of precision
	      double var_new = (pow(omega(0,1),2) + omega(0,0) * omega(1,1))/n;
	      double zscore = omega(0,1)/sqrt(var_new);
	      z_score_precision(i,j) = zscore;
	      z_score_precision(j,i) = z_score_precision(i,j);
	      p_precision(i,j) = 2 * Rf_pnorm5(-abs(zscore), 0.0, 1.0, 1, 0);
	      p_precision(j,i) = p_precision(i,j);

	      // 95% confidence interval of precision
	      double z_95CI = 1.96; // z-score of 95% CI
	      CI_low_precision(i,j) = precision(i,j) - z_95CI * sqrt(var_new);
	      CI_low_precision(j,i) = CI_low_precision(i,j);
	      CI_high_precision(i,j) = precision(i,j) + z_95CI * sqrt(var_new);
	      CI_high_precision(j,i) = CI_high_precision(i,j);

	      // Partial correlation
	      partialCor(i,j) = -omega(0,1)/sqrt(omega(0,0) * omega(1,1));
	      partialCor(j,i) = partialCor(i,j);

	      // P-value of partial correlation
	      var_new = pow((1-pow(partialCor(i,j),2)),2)/n;
	      zscore = partialCor(i,j)/sqrt(var_new);
	      z_score_partialCor(i,j) = zscore;
	      z_score_partialCor(j,i) = z_score_partialCor(i,j);
	      p_partialCor(i,j) = 2 * Rf_pnorm5(-abs(zscore), 0.0, 1.0, 1, 0);
	      p_partialCor(j,i) = p_partialCor(i,j);

	      // 95% confidence interval of partial correlation
	      CI_low_parCor(i,j) = partialCor(i,j) - z_95CI * (1 - pow(partialCor(i,j),2))/sqrt(n);
	      CI_low_parCor(j,i) = CI_low_parCor(i,j);
	      CI_high_parCor(i,j) = partialCor(i,j) + z_95CI * (1 - pow(partialCor(i,j),2))/sqrt(n);
	      CI_high_parCor(j,i) = CI_high_parCor(i,j);
	    }
	  }

	  if(global == true){

	  	    Rcout <<"Perform global inference." << endl;
	  	    Rcout << "Use pre-specified level(s): " << rcpp_level << endl;

	  	    //Take out the upper triangular of test statistic
	  	    size_t t = 0;
	  	    NumericVector upper_stat(p*(p-1)/2);
	  	    for (size_t i = 1; i < p; i++){
	  	    	for (size_t j = 0; j < i; j++){
	  	    		upper_stat[t]=z_score_precision(j,i);
	  	    		t++;
	  	    	}
	  	    }

	  	   //Take out the upper triangular of true graph if available
	  	    NumericMatrix rcpp_true_graph(p,p);
	  	    NumericVector omega_true(p*(p-1)/2);
	  	    t = 0;

	  	    if(true_graph.isNull()){
	  	    	Rcout << "True graph is not available." << endl;
	  	    }
	  	    else{
	  	    	Rcout << "True graph is available." << endl;
	  	       rcpp_true_graph = as<NumericMatrix>(true_graph);
	  		   for (size_t i = 1; i < p; i++){
	  			   for (size_t j = 0; j < i; j++){
	  				   omega_true[t]=rcpp_true_graph(j,i);
	  				   t++;
	  			   }
	  		   }
	  	   }

	  	   // Get threshold
	  	   double upper = 2*sqrt(log(p));
	  	   NumericVector z(2000);
	  	   NumericVector numerator1(2000);
	  	   NumericVector denominator(2000);
	  	   for(double k = 0; k < 2000; k++){
	  		   z[k] = ((k+1)/2000) * upper;
	  		   numerator1[k] = 2 * (1-Rf_pnorm5(z[k], 0.0, 1.0, 1, 0)) *p*(p-1)/2;
	  		   denominator[k] = sum(abs(upper_stat)>=z[k]);
	  		   if(denominator[k] < 1){
	  			   denominator[k] = 1;
	  		   }
	  	   }
	  	   NumericVector FDP = numerator1/denominator;

	  	   NumericVector position(alpha_size);
		   for(size_t i = 0; i < alpha_size; i++){
			   for (double k = 0; k < 2000; k++){
				   if(FDP[k] <= rcpp_level[i]){		// threshold for FDP <= pre-specified level
				   			threshold[i] = z[k];
				   			position[i] = k;
				   			break;
				   }
				   else{
					   if(k == 1999){
						   threshold[i] = z[k];
						   position[i] = k;
					   }
				   }
			   }
		   }


	  	   if(true_graph.isNull()){
			   for(size_t i = 0; i < alpha_size; i++){
			   	   		   FDR[i] = FDP[position[i]];
			   }
		   	   for(size_t k = 0; k < alpha_size; k++){
		   		   NumericMatrix decision(p,p);
		   		   for(size_t i = 0; i < p-1; i++){
		   			   for(size_t j = i+1; j < p; j++){
		   				   if(abs(z_score_precision(i,j)) >= threshold[k]){
		   					   decision(i,j) = 1;
		   				   }
		   				   decision(j,i) = decision(i,j);
		   			   }
		   		   }
		   		   global_decision[k] = decision;
		   	   }
	  	   	}
	  	   	   else{
		   		   for(size_t i = 0; i < alpha_size; i++){
		   			   double FDR_denominator = sum(abs(upper_stat)>= threshold[i]);
		   			   if(FDR_denominator < 1){
		   				   	   FDR_denominator = 1;
		   			   }

		   			   FDR[i] = sum((abs(upper_stat) >= threshold[i]) & (omega_true == 0)) / FDR_denominator;

		   			   double Power_denominator = sum(omega_true!=0);
		   			   if(Power_denominator < 1){
		   				   Power_denominator = 1;
		   			   }
		   			   Power[i] = sum((abs(upper_stat) >= threshold[i]) & (omega_true != 0)) / Power_denominator;
		   		   }
			   	   for(size_t k = 0; k < alpha_size; k++){
			   		   NumericMatrix decision(p,p);
			   		   for(size_t i = 0; i < p-1; i++){
			   			   for(size_t j = i+1; j < p; j++){
			   				   if(abs(z_score_precision(i,j)) >= threshold[k]){
			   					   decision(i,j) = 1;
			   				   }
			   				   decision(j,i) = decision(i,j);
			   			   }
			   		   }
			   		   global_decision[k] = decision;
			   	   }
	  	   	}
	  }

	  if(global == false){
		  return List::create(_["precision"] = precision, _["partialCor"] = partialCor, _["z_score_precision"] = z_score_precision, _["z_score_partialCor"] = z_score_partialCor, _["p_precision"] = p_precision, _["p_partialCor"] = p_partialCor, _["CI_low_precision"] = CI_low_precision, _["CI_high_precision"] = CI_high_precision, _["CI_low_partialCor"] = CI_low_parCor, _["CI_high_partialCor"] = CI_high_parCor);
	  }
	  else{
	  	  if(true_graph.isNull()){
	  		return List::create(_["precision"] = precision, _["partialCor"] = partialCor, _["z_score_precision"] = z_score_precision, _["z_score_partialCor"] = z_score_partialCor, _["p_precision"] = p_precision, _["p_partialCor"] = p_partialCor, _["CI_low_precision"] = CI_low_precision, _["CI_high_precision"] = CI_high_precision, _["CI_low_partialCor"] = CI_low_parCor, _["CI_high_partialCor"] = CI_high_parCor, _["threshold"] = threshold, _["FDR"] = FDR, _["global_decision"] = global_decision);
	  	  }
	  	  else{
	  		return List::create(_["precision"] = precision, _["partialCor"] = partialCor, _["z_score_precision"] = z_score_precision, _["z_score_partialCor"] = z_score_partialCor, _["p_precision"] = p_precision, _["p_partialCor"] = p_partialCor, _["CI_low_precision"] = CI_low_precision, _["CI_high_precision"] = CI_high_precision, _["CI_low_partialCor"] = CI_low_parCor, _["CI_high_partialCor"] = CI_high_parCor, _["threshold"] = threshold, _["FDR"] = FDR, _["power"] = Power, _["global_decision"] = global_decision);
	  	  }

	    }


  }

  return R_NilValue;
}






