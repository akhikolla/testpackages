#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _adpss_exact_est_norm_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _adpss_sample_size_norm_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _adpss_work_test_norm_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_adpss_exact_est_norm_c",   (DL_FUNC) &_adpss_exact_est_norm_c,    9},
    {"_adpss_sample_size_norm_c", (DL_FUNC) &_adpss_sample_size_norm_c,  8},
    {"_adpss_work_test_norm_c",   (DL_FUNC) &_adpss_work_test_norm_c,   21},
    {NULL, NULL, 0}
};

void R_init_adpss(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
