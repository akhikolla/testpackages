#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP wCorr_cont(SEXP, SEXP, SEXP);
extern SEXP wCorr_discord(SEXP);
extern SEXP wCorr_fixxFast(SEXP, SEXP);
extern SEXP wCorr_fscale_cutsFast(SEXP);
extern SEXP wCorr_imapThetaFast(SEXP);
extern SEXP wCorr_imapThetaFast2(SEXP);
extern SEXP wCorr_lnlFast(SEXP, SEXP);
extern SEXP wCorr_mainF(SEXP, SEXP, SEXP, SEXP);
extern SEXP wCorr_mapThetaFast(SEXP);
extern SEXP wCorr_optFcFast(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP wCorr_optFFast(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP wCorr_rcpparma_bothproducts(SEXP);
extern SEXP wCorr_rcpparma_hello_world();
extern SEXP wCorr_rcpparma_innerproduct(SEXP);
extern SEXP wCorr_rcpparma_outerproduct(SEXP);
extern SEXP wCorr_tableFast(SEXP, SEXP, SEXP);
extern SEXP wCorr_theta(SEXP);
extern SEXP wCorr_wrankFast(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"wCorr_cont",                  (DL_FUNC) &wCorr_cont,                  3},
    {"wCorr_discord",               (DL_FUNC) &wCorr_discord,               1},
    {"wCorr_fixxFast",              (DL_FUNC) &wCorr_fixxFast,              2},
    {"wCorr_fscale_cutsFast",       (DL_FUNC) &wCorr_fscale_cutsFast,       1},
    {"wCorr_imapThetaFast",         (DL_FUNC) &wCorr_imapThetaFast,         1},
    {"wCorr_imapThetaFast2",        (DL_FUNC) &wCorr_imapThetaFast2,        1},
    {"wCorr_lnlFast",               (DL_FUNC) &wCorr_lnlFast,               2},
    {"wCorr_mainF",                 (DL_FUNC) &wCorr_mainF,                 4},
    {"wCorr_mapThetaFast",          (DL_FUNC) &wCorr_mapThetaFast,          1},
    {"wCorr_optFcFast",             (DL_FUNC) &wCorr_optFcFast,             6},
    {"wCorr_optFFast",              (DL_FUNC) &wCorr_optFFast,              5},
    {"wCorr_rcpparma_bothproducts", (DL_FUNC) &wCorr_rcpparma_bothproducts, 1},
    {"wCorr_rcpparma_hello_world",  (DL_FUNC) &wCorr_rcpparma_hello_world,  0},
    {"wCorr_rcpparma_innerproduct", (DL_FUNC) &wCorr_rcpparma_innerproduct, 1},
    {"wCorr_rcpparma_outerproduct", (DL_FUNC) &wCorr_rcpparma_outerproduct, 1},
    {"wCorr_tableFast",             (DL_FUNC) &wCorr_tableFast,             3},
    {"wCorr_theta",                 (DL_FUNC) &wCorr_theta,                 1},
    {"wCorr_wrankFast",             (DL_FUNC) &wCorr_wrankFast,             2},
    {NULL, NULL, 0}
};

void R_init_wCorr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
