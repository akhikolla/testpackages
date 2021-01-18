#ifndef BVCUTILITIES_H
#define BVCUTILITIES_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double rtnorm(double mu, double sigma);
double rtnorm0(double mu, double sigma);
double rinvgaussian(double mu, double lambda);
double rinvGauss(double mu, double lambda);
arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma);
arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma, double tol);
#endif
