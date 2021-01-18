#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export(name = "DmaxVec")]]
NumericVector DmaxVec(const NumericVector &Y, const int bint) {
  unsigned int b = bint;
  unsigned int a = Y.size() - 2u * b + 1u;
  double left, right;
  NumericVector ret = NumericVector(a);
  
  for (unsigned int i = 0u; i < a; ++i) {
    left = 0.0;
    right = 0.0;
    
    for (unsigned int j = 0u; j < b; ++j) {
      left += Y[i + j];
      right += Y[i + b + j];
    }
    
    ret[i] = std::abs(left - right) / b;
  }

  return ret;
}
