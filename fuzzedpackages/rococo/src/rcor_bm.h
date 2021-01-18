#ifndef _RCOR_BM_H
#define _RCOR_BM_H

#include <Rcpp.h>
 
RcppExport SEXP rcor_matrix_classical(SEXP vx, SEXP r);
RcppExport SEXP rcor_matrix_linear(SEXP vx, SEXP r);
RcppExport SEXP rcor_matrix_exp(SEXP vx, SEXP r);
RcppExport SEXP rcor_matrix_gauss(SEXP vx, SEXP r);
RcppExport SEXP rcor_matrix_epstol(SEXP vx, SEXP r);

RcppExport SEXP rcor_matrices_classical(SEXP vx, SEXP vy, SEXP r1, SEXP r2);
RcppExport SEXP rcor_matrices_linear(SEXP vx, SEXP vy, SEXP r1, SEXP r2);
RcppExport SEXP rcor_matrices_exp(SEXP vx, SEXP vy, SEXP r1, SEXP r2);
RcppExport SEXP rcor_matrices_gauss(SEXP vx, SEXP vy, SEXP r1, SEXP r2);
RcppExport SEXP rcor_matrices_epstol(SEXP vx, SEXP vy, SEXP r1, SEXP r2);

#endif
