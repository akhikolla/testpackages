#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _RKHSMetaMod_calc_Kv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RKHSMetaMod_grplasso(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RKHSMetaMod_grplasso_q(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RKHSMetaMod_mu_max(SEXP, SEXP);
extern SEXP _RKHSMetaMod_pen_MetMod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RKHSMetaMod_penMetaMod_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RKHSMetaMod_PredErr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RKHSMetaMod_RKHSgrplasso(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RKHSMetaMod_RKHSMetMod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _RKHSMetaMod_RKHSMetMod_qmax(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_RKHSMetaMod_calc_Kv",         (DL_FUNC) &_RKHSMetaMod_calc_Kv,          6},
    {"_RKHSMetaMod_grplasso",        (DL_FUNC) &_RKHSMetaMod_grplasso,         7},
    {"_RKHSMetaMod_grplasso_q",      (DL_FUNC) &_RKHSMetaMod_grplasso_q,       5},
    {"_RKHSMetaMod_mu_max",          (DL_FUNC) &_RKHSMetaMod_mu_max,           2},
    {"_RKHSMetaMod_pen_MetMod",      (DL_FUNC) &_RKHSMetaMod_pen_MetMod,      10},
    {"_RKHSMetaMod_penMetaMod_cpp",  (DL_FUNC) &_RKHSMetaMod_penMetaMod_cpp,  12},
    {"_RKHSMetaMod_PredErr",         (DL_FUNC) &_RKHSMetaMod_PredErr,          8},
    {"_RKHSMetaMod_RKHSgrplasso",    (DL_FUNC) &_RKHSMetaMod_RKHSgrplasso,     5},
    {"_RKHSMetaMod_RKHSMetMod",      (DL_FUNC) &_RKHSMetaMod_RKHSMetMod,       7},
    {"_RKHSMetaMod_RKHSMetMod_qmax", (DL_FUNC) &_RKHSMetaMod_RKHSMetMod_qmax,  9},
    {NULL, NULL, 0}
};

void R_init_RKHSMetaMod(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
