#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP cov_aj(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gen_msm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP los_cp(SEXP, SEXP, SEXP);
extern SEXP los_nocp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"cov_aj",   (DL_FUNC) &cov_aj,   5},
    {"gen_msm",  (DL_FUNC) &gen_msm,  7},
    {"los_cp",   (DL_FUNC) &los_cp,   3},
    {"los_nocp", (DL_FUNC) &los_nocp, 3},
    {NULL, NULL, 0}
};

void R_init_etm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
