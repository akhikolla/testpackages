
#ifndef KB_SHARED_H
#define KB_SHARED_H

#include <Rcpp.h>

// sample from standard uniform distribution

inline double rng_unif() {
  double u;
  // same as in base R
  do {
    u = R::unif_rand();
  } while (u <= 0.0 || u >= 1.0);
  return u;
}

// sample integers from discrete distribution
// cumul_weights is a vector of cumulative weights

inline int sample_int(const Rcpp::NumericVector& cumul_weights) {
  double u = rng_unif();
  int j;
  for (j = 0; j < cumul_weights.length(); j++) {
    if (cumul_weights[j] >= u)
      break;
  }
  return j;
}

#endif
