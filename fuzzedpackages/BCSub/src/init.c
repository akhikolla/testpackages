#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP BCSub_calSim(SEXP);
extern SEXP BCSub_dmvnrm_arma(SEXP, SEXP, SEXP, SEXP);
extern SEXP BCSub_mvrnormArma(SEXP, SEXP, SEXP);
extern SEXP BCSub_myfind(SEXP, SEXP);
extern SEXP BCSub_polyurncpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BCSub_samLamV3Cpp(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"BCSub_calSim",      (DL_FUNC) &BCSub_calSim,      1},
    {"BCSub_dmvnrm_arma", (DL_FUNC) &BCSub_dmvnrm_arma, 4},
    {"BCSub_mvrnormArma", (DL_FUNC) &BCSub_mvrnormArma, 3},
    {"BCSub_myfind",      (DL_FUNC) &BCSub_myfind,      2},
    {"BCSub_polyurncpp",  (DL_FUNC) &BCSub_polyurncpp,  8},
    {"BCSub_samLamV3Cpp", (DL_FUNC) &BCSub_samLamV3Cpp, 5},
    {NULL, NULL, 0}
};

void R_init_BCSub(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
