#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP FSInteract_RIT_1class(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FSInteract_RIT_2class(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"FSInteract_RIT_1class", (DL_FUNC) &FSInteract_RIT_1class,  8},
  {"FSInteract_RIT_2class", (DL_FUNC) &FSInteract_RIT_2class, 11},
  {NULL, NULL, 0}
};

void R_init_FSInteract(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
