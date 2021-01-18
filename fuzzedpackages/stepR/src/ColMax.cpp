#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(name = ".colMax")]]
NumericVector colMax(const NumericMatrix &stat) {
  NumericVector ret = NumericVector(stat.ncol());
  
 for (unsigned int i = 0u; i < static_cast<unsigned int>(stat.ncol()); ++i) {
   ret[i] = max(stat( _, i));
 }
  
  return ret;
}
