#include <Rcpp.h>
#include <cstdint>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double FM_Cpp(const NumericVector& X) {
  int n = X.size();
  double TS_sum = 0;
  int i,j;
  short mult;
  double XimXj, XipXj;

  for(i = 0; i < n; i++) {
    for(j = 0; j <= i; j++) {
      // diagonal elements will be counted once and non diagonal - twice
      mult = ( i == j ? 1 : 2 );
      XimXj = X[i] - X[j];
      XipXj = X[i] + X[j];

      TS_sum += mult * (exp(-0.5*XimXj*XimXj) - exp(-0.5*XipXj*XipXj));
    }
  }

  double TS_value = TS_sum / 2 / n;

  return TS_value;
}

