#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(name = ".convolve")]]
NumericVector convolve(const NumericVector &val, const NumericVector &kern) {
  unsigned int m = kern.size();
  unsigned int n = val.size() - m + 1u;
  NumericVector ret = NumericVector(n);
  
  for (unsigned int i = 0u; i < n; ++i) {
    ret[i] = 0.0;
    for (unsigned int j = 0u, k = i + m - 1u; j < m; ++j, --k) {
      ret[i] += val[k] * kern[j];
    }
  }

  return ret;
}

// [[Rcpp::export(name = ".convolveOversampling")]]
NumericVector convolveOversampling(const NumericVector &val, const NumericVector &kern, unsigned int factor) {
  unsigned int m = kern.size();
  unsigned int n = val.size() / factor - m / factor + 1u;

  NumericVector ret = NumericVector(n);
  
  for (unsigned int i = 0u; i < n; ++i) {
    ret[i] = 0.0;
    for (int j = 0u, k = i * factor + m - 1u; j < m; ++j, --k) {
      ret[i] += val[k] * kern[j];
    }
  }
  
  return ret;
}
