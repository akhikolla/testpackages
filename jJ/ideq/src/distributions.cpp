#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "misc_helpers.h"

using namespace Rcpp;

double ldiwishart(const arma::cube & x, const double df,
                  const arma::cube & scale) {
  const double P = x.n_cols;
  double d = 0;
  
  if (x.n_slices != scale.n_slices)
    throw std::invalid_argument("x and scale must have same number of slices");
  
  for (unsigned int i=0; i<x.n_slices; ++i) {
    d += log(arma::det(scale.slice(i))) * df/2.0;
    d -= log(arma::det(x.slice(i))) * (df+P+1.0)/2.0;
    d -= arma::trace(arma::solve(x.slice(i), scale.slice(i))) / 2.0;
  }
  
  return d;
}

double ldmvnorm(const arma::mat & x, const arma::mat & mu, const arma::mat & Sigma) {
  if (arma::size(x) != arma::size(mu)) {
    throw std::invalid_argument("x and mu must have same size");
  }
  if (!Sigma.is_square()) {
    throw std::invalid_argument("Sigma must be square");
  }
  if (mu.n_elem != Sigma.n_rows) {
    throw std::invalid_argument("Number of elements in x, mu must equal dimension of Sigma");
  }
  
  arma::colvec d = arma::vectorise(x - mu);
  arma::mat tmp = d.t() * arma::solve(Sigma, d);
  return -arma::as_scalar(tmp)/2.0;
}

double rigamma(const double a, const double scl) {
  return (1.0 / R::rgamma(a, 1.0/scl));
}

arma::colvec rmvnorm(const arma::colvec & mean, const arma::mat & Sigma) {
  const int n = mean.n_elem;
  arma::colvec z(n);
  for (int i = 0; i < n; ++i) {
    z.at(i) = R::rnorm(0.0, 1.0);
  }

  // Calculate square root matrix of Sigma
  arma::mat Sigma_sqrt;
  try {
    Sigma_sqrt = arma::sqrtmat_sympd(Sigma);
  }
  catch (std::runtime_error & e)
  { 
    Rcout << "Failed to calculate sqrt(Sigma) using sqrtmat_sympd" << std::endl;
    Rcout << "Forcing symmetry and trying again" << std::endl;
    Sigma_sqrt = forceSqrtMat(Sigma);
  }

  arma::colvec x = mean + Sigma_sqrt * z;
  return x;
}
