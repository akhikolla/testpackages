#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _lowpassFilter_convolve(SEXP, SEXP);
extern SEXP _lowpassFilter_convolveOversampling(SEXP, SEXP, SEXP);
extern SEXP _lowpassFilter_deconvolveJump(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _lowpassFilter_deconvolvePeak(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_lowpassFilter_convolve",             (DL_FUNC) &_lowpassFilter_convolve,              2},
  {"_lowpassFilter_convolveOversampling", (DL_FUNC) &_lowpassFilter_convolveOversampling,  3},
  {"_lowpassFilter_deconvolveJump",       (DL_FUNC) &_lowpassFilter_deconvolveJump,        8},
  {"_lowpassFilter_deconvolvePeak",       (DL_FUNC) &_lowpassFilter_deconvolvePeak,       10},
  {NULL, NULL, 0}
};

void R_init_lowpassFilter(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
