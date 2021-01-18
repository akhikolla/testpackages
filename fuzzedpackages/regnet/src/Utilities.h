#ifndef UTILITIES_h
#define UTILITIES_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double Soft(double z, double lambda);
double validation_LS(const arma::mat& x, const arma::vec& y, const arma::vec& b);
double validation_LAD(const arma::mat& x, const arma::vec& y, const arma::vec& b);
double validation_logit(const arma::mat& x0, const arma::vec& y0, const arma::vec& b);
arma::vec fastLm(const arma::vec & y, const arma::mat & X);

#endif
