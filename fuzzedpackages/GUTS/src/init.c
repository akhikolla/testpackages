#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _GUTS_guts_engine(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GUTS_guts_engine", (DL_FUNC) &_GUTS_guts_engine, 3},
    {NULL, NULL, 0}
};

void R_init_GUTS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
