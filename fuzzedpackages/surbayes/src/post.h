
#ifndef POST_H
#define POST_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;



// Functions
arma::vec fastKronEye_Y ( arma::mat const& Sigma, arma::mat const& Y, int const& n, int const& J );
arma::mat fastKronEye_crossprod (
    arma::mat const& XtX,
    arma::mat const& Sigma,
    arma::vec const& pvec,
    int const& n,
    int const& J
);

#endif
