#include <Rcpp.h>
#include <cstdint>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double BHC1_Cpp(const NumericVector& X, double a) {
  double n = X.size();
  NumericVector aXs = abs(clone(X).sort());
  double TS_sum = 0;
  int i,j;
  double sqdiff;

  for (i = 1; i <= n; i++) {
    TS_sum += (n - 2*i + 1) * (n - 2*i + 1) / a / a;
  }

  for(i = 1; i <= n; i++) {
    for(j = 1; j < i; j++) {
      // non diagonal are counted twice
      sqdiff = (aXs[i-1] - aXs[j-1]) * (aXs[i-1] - aXs[j-1]);
      TS_sum += 2 * (n - 2*i + 1) * (n - 2*j + 1) / (a*a + sqdiff);
    }
  }

  double TS_value = 2 * a / n / (n-1) / (n-1) * TS_sum;

  return TS_value;
}

