#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP estimateParams(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP predict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"estimateParams", (DL_FUNC) &estimateParams, 6},
  {"predict",        (DL_FUNC) &predict,        7},
  {NULL, NULL, 0}
};

void R_init_psgp(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
