#ifndef _COVEM_H
#define _COVEM_H
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
RcppExport SEXP CovEM_Rcpp(SEXP X, SEXP N, SEXP P, SEXP Theta0, SEXP G, SEXP D, 
	SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Miss_group_obs_col, SEXP Miss_group_mis_col,
	SEXP Miss_group_p, SEXP Miss_group_n, SEXP Tol, SEXP Maxiter);
#endif
