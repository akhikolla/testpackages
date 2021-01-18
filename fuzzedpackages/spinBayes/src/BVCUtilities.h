#ifndef BVCUTILITIES_H
#define BVCUTILITIES_H

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec mvrnormChol(const arma::vec& mu, const arma::mat& sigma);
double rinvgaussian(double mu, double lambda);
arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma);
arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma, double tol);
// arma::vec mvrnormEigen(arma::vec& mu, arma::mat& sigma);
#endif
