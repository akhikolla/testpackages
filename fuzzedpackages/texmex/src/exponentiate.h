#ifndef EXPONENTIATE_H
#define EXPONENTIATE_H

#include <Rcpp.h>

namespace {

  template
  <typename T1>
  inline
  Rcpp::NumericVector perhaps_exp(const T1& y, bool log) {
    if (log) {
      return y;
    } else {
      return Rcpp::exp(y);
    }
  }
  
}

#endif
