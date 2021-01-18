#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP clere(SEXP);
extern SEXP pacs(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"clere", (DL_FUNC) &clere, 1},
  {"pacs",  (DL_FUNC) &pacs,  1},
  {NULL, NULL, 0}
};

void R_init_clere(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
