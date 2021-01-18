#include "RcppArmadillo.h"
#include "sample.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @name comire.internal
//' @keywords internal
// [[Rcpp::export(.labelling_c_C)]]
arma::vec labelling_c_C(arma::vec y, arma::vec logpi, arma::vec mu, arma::vec tau){
  
  int n0 = y.n_elem; 
  int H = logpi.n_elem;
  
  arma::vec c(n0);
  // Rcpp::IntegerVector clusters = Rcpp::Range(1,H); 
  arma::vec clusters = arma::regspace(1, H);
  arma::vec lprob(H); 
  arma::vec prob(H);
  
  for(int i=0; i < n0; i++) {
    for(int h =0; h < H; h++) {
      lprob[h] = logpi(h) + R::dnorm(y(i),mu(h),1/sqrt(tau(h)),true);
    }
    lprob = lprob - max(lprob);
    prob  = exp(lprob);
    prob  = prob/sum(prob);
    c[i] = sample(clusters, prob);
  }
  return c;
}
