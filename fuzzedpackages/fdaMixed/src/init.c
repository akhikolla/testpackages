#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP boxcoxTransform(SEXP, SEXP, SEXP, SEXP);
extern SEXP fdaEngine(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fdaTrace(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP marginalANOVA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"boxcoxTransform", (DL_FUNC) &boxcoxTransform,  4},
  {"fdaEngine",       (DL_FUNC) &fdaEngine,       19},
  {"fdaTrace",        (DL_FUNC) &fdaTrace,        10},
  {"marginalANOVA",   (DL_FUNC) &marginalANOVA,   10},
  {NULL, NULL, 0}
};

void R_init_fdaMixed(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
