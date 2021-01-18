#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(name = ".inOrdered")]]
LogicalVector inOrdered(const IntegerVector &x, const IntegerVector &table) {
  LogicalVector ret = LogicalVector(x.size());
  
  int i = 0, j = 0;
  while (i < x.size() && j < table.size()) {
    if (x[i] == table[j]) {
      ret[i] = true;
      ++i;
      ++j;
    } else {
      if (x[i] < table[j]) {
        ++i;
      } else {
        ++j;
      }
    }
  }
  
  return ret;
}
