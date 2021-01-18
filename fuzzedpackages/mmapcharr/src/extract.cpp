/******************************************************************************/

#include <mmapcharr/charsep-acc.h>

using namespace Rcpp;
using std::size_t;

/******************************************************************************/

template <typename T, int RTYPE>
Vector<RTYPE> extractVec(charSepAcc<T, RTYPE> macc, 
                         const IntegerMatrix& elemInd) {
  
  size_t K = elemInd.nrow();
  Vector<RTYPE> res(K);
  
  for (size_t k = 0; k < K; k++)
    res[k] = macc(elemInd(k, 0) - 1, elemInd(k, 1) - 1);
  
  return res;
}

#define EXTRACT_VEC(CTYPE, RTYPE) {                                            \
return extractVec<CTYPE, RTYPE>(charSepAcc<CTYPE, RTYPE>(xpMat, e["code"]),    \
                                elemInd);                                      \
}

// Dispatch function for extractVec
// [[Rcpp::export]]
RObject extractVec(Environment e,
                   const IntegerMatrix& elemInd) {
  
  XPtr<charSep> xpMat = e["address"];
  
  switch(TYPEOF(e["code"])) {
  case INTSXP:
    EXTRACT_VEC(int,           INTSXP)
  case REALSXP:
    EXTRACT_VEC(double,        REALSXP)
  case RAWSXP:
    EXTRACT_VEC(unsigned char, RAWSXP)
  case STRSXP:
    EXTRACT_VEC(String,        STRSXP)
  case LGLSXP:
    EXTRACT_VEC(int,           LGLSXP)
  default:
    throw Rcpp::exception("Type not supported.");
  }
}

/******************************************************************************/

template <typename T, int RTYPE>
Matrix<RTYPE> extractMat(charSepAcc<T, RTYPE> macc,
                         const IntegerVector& rowInd,
                         const IntegerVector& colInd) {
  
  size_t n = rowInd.size();
  size_t m = colInd.size();
  
  IntegerVector rows = rowInd - 1;
  IntegerVector cols = colInd - 1;
  
  Matrix<RTYPE> res(n, m);
  
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < m; j++)
      res(i, j) = macc(rows[i], cols[j]);
  
  return res;
}

#define EXTRACT_MAT(CTYPE, RTYPE) {                                            \
return extractMat<CTYPE, RTYPE>(charSepAcc<CTYPE, RTYPE>(xpMat, e["code"]),    \
                                rowInd, colInd);                               \
}

// Dispatch function for extractMat
// [[Rcpp::export]]
RObject extractMat(Environment e,
                   const IntegerVector& rowInd,
                   const IntegerVector& colInd) {
  
  XPtr<charSep> xpMat = e["address"];
  
  switch(TYPEOF(e["code"])) {
  case INTSXP:
    EXTRACT_MAT(int,           INTSXP)
  case REALSXP:
    EXTRACT_MAT(double,        REALSXP)
  case RAWSXP:
    EXTRACT_MAT(unsigned char, RAWSXP)
  case STRSXP:
    EXTRACT_MAT(String,        STRSXP)
  case LGLSXP:
    EXTRACT_MAT(int,           LGLSXP)
  default:
    throw Rcpp::exception("Type not supported.");
  }
}

/******************************************************************************/
