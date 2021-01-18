#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _gRain_propagateLS__(SEXP, SEXP);
extern SEXP _gRain_sparse_setXtf1(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_gRain_propagateLS__",  (DL_FUNC) &_gRain_propagateLS__,  2},
    {"_gRain_sparse_setXtf1", (DL_FUNC) &_gRain_sparse_setXtf1, 2},
    {NULL, NULL, 0}
};

void R_init_gRain(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
