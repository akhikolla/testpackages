#include <Rcpp.h>
#include <cstdint>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double SGN_Cpp(const NumericVector& X) {
  int n = X.size();
  return sqrt((double)n) * (mean((X) > 0) - 0.5);
}
