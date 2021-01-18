#include <Rcpp.h>
#include <cstdint>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
double BHI_Cpp(const NumericVector & X) {
  int n = X.size();
  NumericVector aXs = abs(clone(X).sort());

  IntegerVector ranks = match(aXs, clone(aXs).sort());

  double T2 = 0, T3 = 0;

  // We're using 1-indexing
  for (int j = 1; j <= n; j++) {
    T2 += (n - ranks[j - 1]) * Rf_choose(n - j, 1);

    T3 += (n - ranks[j - 1]) *
      Rf_choose(j - 1, 1);
  }

  return sqrt((double)n) * (T2 - T3) / 2 / (n * Rf_choose(n, 2));
}

