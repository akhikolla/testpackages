#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP adkhi2partial(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP computeProba(SEXP, SEXP, SEXP, SEXP);
extern SEXP freqMultiR(SEXP, SEXP);
extern SEXP kullback(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP loglikelihood(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP semR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simulISRR(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"adkhi2partial", (DL_FUNC) &adkhi2partial,  5},
  {"computeProba",  (DL_FUNC) &computeProba,   4},
  {"freqMultiR",    (DL_FUNC) &freqMultiR,     2},
  {"kullback",      (DL_FUNC) &kullback,       7},
  {"loglikelihood", (DL_FUNC) &loglikelihood,  9},
  {"semR",          (DL_FUNC) &semR,          12},
  {"simulISRR",     (DL_FUNC) &simulISRR,      4},
  {NULL, NULL, 0}
};

void R_init_Rankcluster(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
