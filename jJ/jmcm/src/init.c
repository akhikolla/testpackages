#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP MCD__new(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ACD__new(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP HPC__new(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP get_m(SEXP, SEXP);
extern SEXP get_Y(SEXP, SEXP);
extern SEXP get_X(SEXP, SEXP);
extern SEXP get_Z(SEXP, SEXP);
extern SEXP get_W(SEXP, SEXP);
extern SEXP get_D(SEXP, SEXP, SEXP);
extern SEXP get_T(SEXP, SEXP, SEXP);
extern SEXP get_mu(SEXP, SEXP, SEXP);
extern SEXP get_Sigma(SEXP, SEXP, SEXP);
extern SEXP n2loglik(SEXP, SEXP);
extern SEXP grad(SEXP, SEXP);
extern SEXP hess(SEXP, SEXP);
extern SEXP _jmcm_mcd_estimation(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _jmcm_acd_estimation(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _jmcm_hpc_estimation(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"MCD__new",             (DL_FUNC) &MCD__new,              5},
    {"ACD__new",             (DL_FUNC) &ACD__new,              5},
    {"HPC__new",             (DL_FUNC) &HPC__new,              5},
    {"get_m",                (DL_FUNC) &get_m,                 2},
    {"get_Y",                (DL_FUNC) &get_Y,                 2},
    {"get_X",                (DL_FUNC) &get_X,                 2},
    {"get_Z",                (DL_FUNC) &get_Z,                 2},
    {"get_W",                (DL_FUNC) &get_W,                 2},
    {"get_D",                (DL_FUNC) &get_D,                 3},
    {"get_T",                (DL_FUNC) &get_T,                 3},
    {"get_mu",               (DL_FUNC) &get_mu,                3},
    {"get_Sigma",            (DL_FUNC) &get_Sigma,             3},
    {"n2loglik",             (DL_FUNC) &n2loglik,              2},
    {"grad",                 (DL_FUNC) &grad,                  2},
    {"hess",                 (DL_FUNC) &hess,                  2},
    {"_jmcm_mcd_estimation", (DL_FUNC) &_jmcm_mcd_estimation, 12},
    {"_jmcm_acd_estimation", (DL_FUNC) &_jmcm_acd_estimation, 12},
    {"_jmcm_hpc_estimation", (DL_FUNC) &_jmcm_hpc_estimation, 12},
    {NULL, NULL, 0}
};

void R_init_jmcm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
