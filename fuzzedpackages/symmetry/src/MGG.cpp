#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MGG_Cpp(const NumericVector& X) {
  return sqrt((double)X.size()) * (mean(X)-median(X)) /
    (sqrt(M_PI/2) * mean(abs(X - median(X))) * sqrt(0.5707963));
}

