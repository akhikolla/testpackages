#include "R.h"
#include "Rinternals.h"
#include "stdlib.h" // for NULL
#include "R_ext/Rdynload.h"

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP basadFunctionG(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP );
extern SEXP basadFuncScale( SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"basadFunctionG", (DL_FUNC) &basadFunctionG, 14},
    {"basadFuncScale", (DL_FUNC) &basadFuncScale, 16},
    {NULL, NULL, 0}
};

void R_init_basad(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
