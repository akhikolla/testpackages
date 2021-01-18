#ifndef _RCOR_H
#define _RCOR_H

#include <Rcpp.h>

using namespace Rcpp;

void get_sums (NumericMatrix mat_x, NumericMatrix mat_y, IntegerVector perm,
	       double (*tnorm_fp)(double, double), double *sum, double *sum_t);

RcppExport SEXP rcor (SEXP matx, SEXP maty, SEXP tnorm);

#endif
