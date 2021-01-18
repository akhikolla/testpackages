#ifndef ROBUST_h
#define ROBUST_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec LadNet(const arma::mat& x, const arma::vec& y, double lam1, double lam2, arma::vec b, double r, const arma::mat& a, int n, int p);
arma::vec LadMCP(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, double r, int n, int p);
arma::vec LadLasso(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, int n, int p);

#endif
