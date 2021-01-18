#include <Rcpp.h>
#include <cstdint>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double BH2_Cpp(const NumericVector& X) {
  double n = X.size();

  double TS_sum = 0;

  // Get the equivalent of order()
  NumericVector aX = abs(X);
  NumericVector sorted = clone(aX).sort();
  IntegerVector order = match(sorted, aX) - 1;

  NumericVector S(n);

  for (int i = 0; i < n; i++) {
    S[i] = (X[order[i]] > 0)? 1 : -1;
    if (X[order[i]] == 0)
      S[i] = 0;
  }

  NumericVector inv_cumsum(n);
  double sum = 0;
  for (int i = n - 1; i >= 0; i--) {
    sum += S[i];
    inv_cumsum[i] = sum;
  }

  for (int i = 0; i < n-1; i++) {
    double inner_sum = inv_cumsum[i+1];
    int k = i+1;
    TS_sum += inner_sum * inner_sum * k * k / n / n;
  }

  double TS_value = TS_sum / (n-1) / (n-1);

  return TS_value;
}

