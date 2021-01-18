#ifndef ScoreU_h
#define ScoreU_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List ScoreU(int n, arma::vec& k, arma::vec& y, arma::mat& x, int p, arma::vec& beta, arma::cube& Rhat);

#endif
