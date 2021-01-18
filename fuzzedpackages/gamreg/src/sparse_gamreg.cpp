//'  @useDynLib(gamreg, .registration = TRUE)
//'  @importFrom Rcpp evalCpp
// [[Rcpp::depends(RcppArmadillo)]]

// RegisteringDynamic Symbols

#include <RcppArmadillo.h>
#include <R.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

//
// [[Rcpp::export]]
List gam_reg( arma::mat X, arma::mat Y, arma::mat beta, double beta0, double sigma, double lambda, double gam, Rcpp::Function f , int inter, double regul_alp  ){

  if(all(beta.col(0)==0)){
    List res;
    res["beta0"] = beta0;
    res["beta"] = beta;
    res["sigma"] =sigma;
    return res;

  }


  const int N = X.n_rows;
  //const int p = X.n_cols;

  const arma::mat tmp = arma::ones(N,1);
  arma::mat tmp1 = inter*beta0*tmp;
  arma::mat tmp2 = X*beta;
  arma::mat tmp3 = tmp1+tmp2;




  for(int i=0; i < 5000 ; i++){
    double beta0_tmp = beta0*inter;
    mat beta_tmp = beta;
    double sigma_tmp =sigma;


    mat alpha = exp(-gam*pow(Y-tmp3,2)/(2*pow(sigma,2)));
    alpha = alpha/accu(alpha);




    beta0 = as_scalar(dot(alpha,Y-tmp2))*inter;

    tmp1 = beta0*tmp;


    mat Y_lasso = diagmat(sqrt(alpha))*(Y-tmp1)/sigma;
    mat  X_lasso =  diagmat(sqrt(alpha))*X/sigma;


    List res = f(X_lasso,Y_lasso, _["intercept"]=0 , _["standardize"]=1,  _["alpha"]=regul_alp, _["family"]="gaussian" , _["lambda"]=lambda);
    sp_mat t =res["beta"];
    mat z(t);
    beta = z;

    tmp2 = X*beta;
    tmp3 = tmp1+tmp2;


    sigma = sqrt((1+gam)*as_scalar(dot(alpha,pow(Y-tmp3,2))));



    if(  all(beta.col(0)==0) ||  (-gam/(1+gam))*log(sigma_tmp) -(-gam/(1+gam))*log(sigma)  + (-1/gam)*log( accu( exp(-gam*pow(Y-beta0_tmp*tmp-X*beta_tmp,2)/(2*pow(sigma_tmp,2)) )* pow( 2*arma::datum::pi*pow(sigma_tmp,2) , -gam/2)) )  + lambda*accu(abs(beta_tmp))- lambda*accu(abs(beta)) - (-1/gam)*log( accu( exp(-gam*pow(Y-tmp3,2)/(2*pow(sigma,2)) )* pow( 2*arma::datum::pi*pow(sigma,2) , -gam/2)) )  <= pow(10,-5) ){
      break;
    }



  }

  List res;
  res["beta0"] = beta0;
  res["beta"] = beta;
  res["sigma"] =sigma;
  return res;

}


void R_init_markovchain(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
