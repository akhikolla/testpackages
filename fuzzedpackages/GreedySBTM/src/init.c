#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* TO GENERATE THIS FILE I HAVE USED:
   tools::package_native_routine_registration_skeleton(".")
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GreedySBTM_cpp_CollapseLabels(SEXP);
extern SEXP GreedySBTM_cpp_GreedyICL(SEXP, SEXP, SEXP, SEXP);
extern SEXP GreedySBTM_cpp_GreedyMerge(SEXP, SEXP, SEXP);
extern SEXP GreedySBTM_cpp_ICLExact(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"GreedySBTM_cpp_CollapseLabels", (DL_FUNC) &GreedySBTM_cpp_CollapseLabels, 1},
    {"GreedySBTM_cpp_GreedyICL",      (DL_FUNC) &GreedySBTM_cpp_GreedyICL,      4},
    {"GreedySBTM_cpp_GreedyMerge",    (DL_FUNC) &GreedySBTM_cpp_GreedyMerge,    3},
    {"GreedySBTM_cpp_ICLExact",       (DL_FUNC) &GreedySBTM_cpp_ICLExact,       3},
    {NULL, NULL, 0}
};

void R_init_GreedySBTM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
