#include <Rcpp.h>
#include <cstdint>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double MOI_Cpp(const NumericVector & X, double k) {
  int n = X.size();
  NumericVector aXs = abs(clone(X).sort());

  IntegerVector ranks = match(aXs, clone(aXs).sort());

  double T2 = 0, T3 = 0;

  // We're using 1-indexing
  for (int j = k; j <= n - k + 1; j++) {
    T2 += (n - ranks[j - 1]) *
      Rf_choose(n - j, k) *
      Rf_choose(j - 1, k - 1);

    T3 += (n - ranks[j - 1]) *
      Rf_choose(n - j, k - 1) *
      Rf_choose(j - 1, k);
  }

  return sqrt((double)n) * (T2 - T3) / (n * Rf_choose(n, 2 * k));
}

