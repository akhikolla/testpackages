#include <Rcpp.h>
#include <algorithm> 

using namespace Rcpp;

// [[Rcpp::export(name = ".kernelSmoothing")]]
NumericVector kernelSmoothing(const NumericVector &Y, const NumericVector &K) {
  int n = Y.size();
  int b = (K.size() - 1u) / 2u;
  double est, sumK;
  NumericVector ret = NumericVector(n);
  
  for (int i = 0u; i < n; ++i) {
    est = 0.0;
    sumK = 0.0;
    
    for (unsigned int j = std::max(i - b, 0), k = std::max(0, b - i); j <= std::min(i + b, n - 1); ++j, ++k) {
      est += Y[j] * K[k];
      sumK += K[k];
    }
    ret[i] = est / sumK;
  }
  
  return ret;
}
