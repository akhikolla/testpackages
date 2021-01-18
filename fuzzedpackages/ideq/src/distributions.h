#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double ldiwishart(const arma::cube & x, const double df,
                  const arma::cube & scale);

double ldmvnorm(const arma::mat & x, const arma::mat & mu, const arma::mat & Sigma);

double rigamma(const double a, const double scl);

arma::colvec rmvnorm(const arma::colvec & mean, const arma::mat & Sigma);

#endif
