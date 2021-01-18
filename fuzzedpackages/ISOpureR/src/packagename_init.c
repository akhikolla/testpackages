#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ISOpureR_rcppeigen_max_over_columns_or_rows(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ISOpureR_rcppeigen_max_over_columns_or_rows", (DL_FUNC) &_ISOpureR_rcppeigen_max_over_columns_or_rows, 2},
    {NULL, NULL, 0}
};

void R_init_ISOpureR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
