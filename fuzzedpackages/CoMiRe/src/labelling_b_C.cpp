#include "RcppArmadillo.h"
#include "sample.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @name comire.internal
//' @keywords internal
// [[Rcpp::export(.labelling_b_C)]]
arma::vec labelling_b_C(arma::vec w, arma::mat phi, arma::vec f0i, arma::vec f1i){
  int n = f0i.n_elem;
  int J = phi.n_cols;

  arma::vec b(n);
  // Rcpp::IntegerVector clusters = Rcpp::Range(1,J); 
  arma::vec clusters = arma::regspace(1, J);
  arma::vec prob(J);

  for(int i=0; i < n; i++) {
    for(int j = 0; j < J; j++) {
      prob(j) = w(j)*((1 - phi(i,j))*f0i(i) + phi(i,j)*f1i(i));
      if (prob(j) < 0) {
        prob(j) = 0;
      }
    }
    prob  = prob/sum(prob);
    b(i) = sample(clusters, prob);
  }
  return b;
}
