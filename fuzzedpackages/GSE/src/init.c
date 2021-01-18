#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP CovEM_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP emve_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fast_partial_mahalanobis(SEXP, SEXP, SEXP, SEXP);
extern SEXP GSE_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GSE_Rocke(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"CovEM_Rcpp",               (DL_FUNC) &CovEM_Rcpp,               14},
    {"emve_Rcpp",                (DL_FUNC) &emve_Rcpp,                21},
    {"fast_partial_mahalanobis", (DL_FUNC) &fast_partial_mahalanobis,  4},
    {"GSE_Rcpp",                      (DL_FUNC) &GSE_Rcpp,            14},
    {"GSE_Rocke",                (DL_FUNC) &GSE_Rocke,                15},
    {NULL, NULL, 0}
};

void R_init_GSE(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
