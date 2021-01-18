#ifndef _EMVE_H
#define _EMVE_H
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
RcppExport SEXP emve_Rcpp(SEXP X, SEXP X_nonmiss, SEXP Pu, SEXP N, SEXP P, SEXP Theta0, SEXP GG, SEXP D, SEXP X_miss_group_match,
    SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Miss_group_obs_col, SEXP Miss_group_mis_col,
    SEXP Miss_group_p, SEXP Miss_group_n, SEXP NResample, SEXP NSubsampleSize, SEXP Subsample_ID, SEXP CC, SEXP CK, SEXP Maxits);
RcppExport SEXP fast_partial_mahalanobis(SEXP X_mu_diff, SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts);
#endif
