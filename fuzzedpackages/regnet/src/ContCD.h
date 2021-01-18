#ifndef CONTCD_h
#define CONTCD_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec ContNet(const arma::mat& x, const arma::vec& y, double lam1, double lam2, arma::vec b, double r, const arma::mat& a, int n, int p);
arma::vec ContMCP(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, double r, int n, int p);
arma::vec ContLasso(const arma::mat& x, const arma::vec& y, double lam1, arma::vec b, int n, int p);

#endif
