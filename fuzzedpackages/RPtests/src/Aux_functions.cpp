#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector colMax(NumericMatrix mat, bool na_rm=false) {
  NumericVector out = mat(0, _);
  if (na_rm) {
    for (int j=0; j<out.size(); j++) {
      if (NumericVector::is_na(out[j])) out[j] = R_NegInf;
    }
    for (int j=0; j<mat.ncol(); j++) {
      for (int i=1; i<mat.nrow(); i++) {
        if (!NumericVector::is_na(mat(i, j))) {
          if (mat(i, j) > out[j]) out[j] = mat(i, j);
        }
      }
    }
  }
  else {
    for (int j=0; j<mat.ncol(); j++) {
      for (int i=1; i<mat.nrow(); i++) {
        if (mat(i, j) > out[j]) out[j] = mat(i, j);
      }
    }
  }
  return out;
}
