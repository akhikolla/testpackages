#ifndef _GSE_Rocke_H
#define _GSE_Rocke_H
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
RcppExport SEXP GSE_Rocke(SEXP X, SEXP N, SEXP P, SEXP Mu0, SEXP S0, SEXP Tol, SEXP Maxiter, SEXP Tol_scale, SEXP Miter_scale, SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Tuning_const_group, SEXP Gamma_tune_group, SEXP Print_step, SEXP Bdp);
#endif
