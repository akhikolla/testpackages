
#include <Rcpp.h>
#include "kernels.h"
#include "shared.h"

using Rcpp::NumericVector;


// [[Rcpp::export]]
NumericVector cpp_ruvk(
    const int& n,
    const NumericVector& y,
    const double& bandwidth,
    const NumericVector& weights,
    const std::string& kernel = "gaussian",
    const bool& shrinked = false
  ) {

  if (y.length() < 1) {
    Rcpp::warning("NAs produced");
    NumericVector out(n, NA_REAL);
    out.attr("boot_index") = NumericVector(n, NA_REAL);
    return out;
  }

  double (*rng_kern)();

  if (kernel == "rectangular") {
    rng_kern = rng_rect;
  } else if (kernel == "triangular") {
    rng_kern = rng_triang;
  } else if (kernel == "biweight") {
    rng_kern = rng_biweight;
  } else if (kernel == "cosine") {
    rng_kern = rng_cosine;
  } else if (kernel == "optcosine") {
    rng_kern = rng_optcos;
  } else if (kernel == "epanechnikov") {
    rng_kern = rng_epan;
  } else {
    rng_kern = R::norm_rand;
  }

  const int k = y.length();
  NumericVector samp(n);
  NumericVector c_weights(k);
  std::vector<int> idx(n);

  if (!R_FINITE(bandwidth))
    Rcpp::stop("inappropriate value of bandwidth");

  if (bandwidth < 0.0)
    Rcpp::stop("bandwidth needs to be non-negative");

  if (Rcpp::is_true(Rcpp::any(weights < 0.0)))
    Rcpp::stop("weights need to be non-negative");

  if (Rcpp::is_false(Rcpp::all(Rcpp::is_finite(weights))))
    Rcpp::stop("inappropriate values of weights");

  if (weights.length() == 1) {
    c_weights.fill( 1.0/static_cast<double>(k) );
  } else {
    if (weights.length() != k)
      Rcpp::stop("dimmensions of weights and y do not match");
    c_weights = weights;
  }

  for (int i = 1; i < k; i++)
    c_weights[i] += c_weights[i-1];
  c_weights = c_weights / c_weights[k-1];

  if (k == 1) {

    for (int i = 0; i < n; i++) {
      samp[i] = y[0] + rng_kern() * bandwidth;
    }

  } else if (!shrinked) {

    int j;
    for (int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      samp[i] = y[j] + rng_kern() * bandwidth;
    }

  } else {

    const double my = Rcpp::mean(y);
    const double sy = Rcpp::var(y);
    const double c = std::sqrt(1.0 + (bandwidth*bandwidth)/sy);

    int j;
    for (int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      samp[i] = my + (y[j] - my + rng_kern() * bandwidth) / c;
    }

  }

  for (int i = (k-1); i > 0; i--)
    c_weights[i] -= c_weights[i-1];

  samp.attr("boot_index") = idx;
  return samp;

}

