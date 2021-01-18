#include <Rcpp.h>
#include <cstdint>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double RW_Cpp(const NumericVector& X) {
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

  for (int i = 0; i < n; i++) {
    double inner_sum = (i != n - 1) ? inv_cumsum[i+1] : 0;
    TS_sum += (S[i] + 2 * inner_sum) * (S[i] + 2 * inner_sum);
  }

  double TS_value = TS_sum / 4 / n / n;

  return TS_value;
}

