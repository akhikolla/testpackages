#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double B1_Cpp(const NumericVector& X) {
  NumericVector Xc = X - mean(X);
  double std = sd(X);
  return mean(Xc*Xc*Xc) / (std*std*std);
}

