#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @name comire.internal
//' @keywords internal
// [[Rcpp::export(.mixdensity_C)]]
arma::vec mixdensity_C(arma::vec y, arma::vec pi, arma::vec mu, arma::vec tau){
  
  int n = y.n_elem; 
  int H = pi.n_elem;
  
  arma::vec f0i(n);
  Rcpp::IntegerVector clusters = Rcpp::Range(1,H); 
  arma::vec prob(H);
  
  for(int i=0; i < n; i++) {
    double prob = 0;
    for(int h =0; h < H; h++) {
      prob = prob + pi(h)*R::dnorm(y(i),mu(h),1/sqrt(tau(h)),false);
    }
    f0i(i) = prob;
  }
  return f0i;
}
