/******************************************************************************/

#include <mmapcharr/charsep-acc.h>
using namespace Rcpp;

/******************************************************************************/

// [[Rcpp::export]]
int nline_cpp(SEXP obj) {
  
  XPtr<charSep> xpMat(obj);
  const unsigned char* mat = xpMat->matrix();
  
  const unsigned char backslash_n = 10;
  
  return std::count(mat, mat + xpMat->nrow(), backslash_n);
}

/******************************************************************************/
