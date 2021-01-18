#include <Rcpp.h>
#include <algorithm> 

using namespace Rcpp;

// [[Rcpp::export(name = ".CVonesided")]]
double CVonesided(const NumericVector &Y, const NumericVector &K) {
  unsigned int n = Y.size();
  unsigned int b = K.size();
  double cv = 0.0, est, sumK;
  
  for (unsigned int i = 1u; i < n; ++i) {
    est = 0.0;
    sumK = 0.0;

    for (unsigned int j = i - 1u, k = 0u; k < std::min(b, i); --j, ++k) {
      est += Y[j] * K[k];
      sumK += K[k];
    }
    cv += (est / sumK - Y[i]) * (est / sumK - Y[i]);
  }

  return cv;
}

// [[Rcpp::export(name = ".CVtwosided")]]
double CVtowsided(const NumericVector &Y, const NumericVector &K) {
  unsigned int n = Y.size();
  unsigned int b = K.size();
  double cv = 0.0, est, sumK;
  
  for (unsigned int i = 0u; i < n; ++i) {
    est = 0.0;
    sumK = 0.0;
    
    for (unsigned int j = i - 1u, k = 0u; k < std::min(b, i); --j, ++k) {
      est += Y[j] * K[k];
      sumK += K[k];
    }
    
    for (unsigned int j = i + 1u, k = 0u; k < std::min(b, n - i - 1u); ++j, ++k) {
      est += Y[j] * K[k];
      sumK += K[k];
    }
    cv += (est / sumK - Y[i]) * (est / sumK - Y[i]);
  }
  
  return cv;
}
