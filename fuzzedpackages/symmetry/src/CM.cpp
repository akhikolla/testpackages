#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double CM_Cpp(const NumericVector& X) {
  return sqrt((double)X.size()) * (mean(X)-median(X)) / (sd(X) * sqrt(0.5707963));
}

