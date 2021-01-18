// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::vec OLS_cpp (const arma::mat& X, const arma::vec& y) {
  return((X.t() * X ).i() * X.t() * y) ;
}                  

// [[Rcpp::export]]
arma::vec OLS2_cpp (const arma::mat& X, const arma::vec& y) {
  return(arma::solve(X, y)) ;
}                  
