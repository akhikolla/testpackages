
#include <Rcpp.h>
#include "kernels.h"
#include "shared.h"

using Rcpp::NumericVector;
using Rcpp::NumericMatrix;


// [[Rcpp::export]]
NumericMatrix cpp_rmvk(
    const int& n,
    const NumericMatrix& y,
    const NumericVector& bandwidth,
    const NumericVector& weights,
    const std::string& kernel = "gaussian",
    const bool& shrinked = false
  ) {

  if (y.nrow() < 1 || y.ncol() < 1) {
    Rcpp::warning("NAs produced");
    NumericMatrix out(n, y.ncol());
    if (y.ncol() < 1)
      return out;
    std::fill(out.begin(), out.end(), NA_REAL);
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

  const int m = y.ncol();
  const int k = y.nrow();
  NumericMatrix samp(n, m);
  NumericVector c_weights(k);
  std::vector<int> idx(n);

  if (bandwidth.length() != m)
    Rcpp::stop("dimmensions of y and bandwidth do not match");

  if (Rcpp::is_false(Rcpp::all(Rcpp::is_finite(bandwidth))))
    Rcpp::stop("inappropriate values of bandwidth");

  if (Rcpp::is_true(Rcpp::any(bandwidth < 0.0)))
    Rcpp::stop("bandwidth needs to be non-negative");

  if (!Rcpp::is_true(Rcpp::any(Rcpp::is_finite(weights))))
    Rcpp::stop("inappropriate values of weights");

  if (Rcpp::is_true(Rcpp::any(weights < 0.0)))
    Rcpp::stop("weights need to be non-negative");

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
      for (int l = 0; l < m; l++)
        samp(i, l) = y(0, l) + rng_kern() * bandwidth[l];
    }

  } else if (!shrinked) {

    int j;
    for (int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      for (int l = 0; l < m; l++)
        samp(i, l) = y(j, l) + rng_kern() * bandwidth[l];
    }

  } else {

    NumericVector my(m);
    for (int i = 0; i < m; i++)
      my[i] = Rcpp::mean(y.column(i));

    NumericVector sy(m);
    for (int i = 0; i < m; i++)
      sy[i] = Rcpp::var(y.column(i));

    NumericVector bw_sq = bandwidth * bandwidth;
    const NumericVector c = Rcpp::sqrt(1.0 + bw_sq/sy);

    int j;
    for (int i = 0; i < n; i++) {
      j = sample_int(c_weights);
      idx[i] = j + 1;
      for (int l = 0; l < m; l++)
        samp(i, l) = my[l] + (y(j, l) - my[l] + rng_kern() * bandwidth[l]) / c[l];
    }

  }

  for (int i = (k-1); i > 0; i--)
    c_weights[i] -= c_weights[i-1];

  samp.attr("boot_index") = idx;
  return samp;

}

