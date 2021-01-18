#ifndef MISC_HELPERS_H
#define MISC_HELPERS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

void makeSymmetric(arma::mat & X);

arma::mat forceInv(arma::mat X);

arma::mat forceSqrtMat(arma::mat X);

#endif
