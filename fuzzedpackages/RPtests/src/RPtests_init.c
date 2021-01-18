#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP RPtests_colMax(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"RPtests_colMax", (DL_FUNC) &RPtests_colMax, 2},
  {NULL, NULL, 0}
};

void R_init_RPtests(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
