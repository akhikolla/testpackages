//=================================
// include guard
#pragma once
#include <RcppArmadillo.h>

const double log2pi = std::log(2.0 * M_PI);

// Multiple outputs
arma::colvec dmvnormMultiple(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd = false);
arma::mat rmvnormMultiple(int n, arma::rowvec mean, arma::mat sigma);

// Single output
double dmvnormSingle(arma::colvec x,  arma::colvec mean,  arma::mat sigma, bool logd = false);
arma::colvec rmvnormSingle(arma::colvec mean, arma::mat sigma);

bool isPositiveDefinite(arma::mat matrix, double tol = 0.00005); 
