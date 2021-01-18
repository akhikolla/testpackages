#ifndef _dils_RSS_CPP_MATRIX_H
#define _dils_RSS_CPP_MATRIX_H

#include <Rcpp.h>
#include <time.h>

RcppExport SEXP rss_cell(SEXP xadj, SEXP vin, SEXP vout, SEXP radius, SEXP directed);
RcppExport SEXP rss_cpp_matrix(SEXP xadj, SEXP radius, SEXP directed);

#endif
