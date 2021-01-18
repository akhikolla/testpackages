#include <Rcpp.h>
#include <iostream>
#include <ctime>

using namespace Rcpp;

template<int RTYPE>
Vector<RTYPE> colsums_sq(Matrix<RTYPE> x) {
  int n = x.nrow(), m = x.ncol();
  Vector<RTYPE> r(m);
  for(int j = 0; j < m ; j++) {
    for(int i = 0; i < n; i++) r(j) += x(i,j)*x(i,j);
  }
  return r;
}

template<int RTYPE>
Vector<RTYPE> colsums_cub(Matrix<RTYPE> x) {
  int n = x.nrow(), m = x.ncol();
  Vector<RTYPE> r(m);
  for(int j = 0; j < m ; j++) {
    for(int i = 0; i < n; i++) r(j) += x(i,j)*x(i,j)*x(i,j);
  }
  return r;
}

SEXP colSumsSq(SEXP x) {
  if(TYPEOF(x) == INTSXP)
    return colsums_sq<INTSXP>(x);
  if(TYPEOF(x) == REALSXP)
    return colsums_sq<REALSXP>(x);
  stop("Not an integer or a numeric matrix");
}

SEXP colSumsCub(SEXP x) {
  if(TYPEOF(x) == INTSXP)
    return colsums_cub<INTSXP>(x);
  if(TYPEOF(x) == REALSXP)
    return colsums_cub<REALSXP>(x);
  stop("Not an integer or a numeric matrix");
}

RcppExport SEXP oz_colsums_sq(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(colSumsSq(x));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP oz_colsums_cub(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(colSumsCub(x));
    return rcpp_result_gen;
END_RCPP
}

