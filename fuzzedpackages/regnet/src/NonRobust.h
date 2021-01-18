#ifndef NONROBUST_h
#define NONROBUST_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec LSNet(const arma::mat& x, const arma::vec& y, double lam1, double lam2, arma::vec b, double r, const arma::mat& a, int n, int p);
arma::vec LSMCP(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, double r, int n, int p);
arma::vec LSLasso(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, int n, int p);

#endif
