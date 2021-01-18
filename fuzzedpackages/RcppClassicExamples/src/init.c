#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP classicRcppDateExample(SEXP, SEXP);
extern SEXP classicRcppMatrixExample(SEXP);
extern SEXP classicRcppParamsExample(SEXP);
extern SEXP classicRcppStringVectorExample(SEXP);
extern SEXP classicRcppVectorExample(SEXP);
extern SEXP Rcpp_Example(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"classicRcppDateExample",         (DL_FUNC) &classicRcppDateExample,         2},
    {"classicRcppMatrixExample",       (DL_FUNC) &classicRcppMatrixExample,       1},
    {"classicRcppParamsExample",       (DL_FUNC) &classicRcppParamsExample,       1},
    {"classicRcppStringVectorExample", (DL_FUNC) &classicRcppStringVectorExample, 1},
    {"classicRcppVectorExample",       (DL_FUNC) &classicRcppVectorExample,       1},
    {"Rcpp_Example",                   (DL_FUNC) &Rcpp_Example,                   9},
    {NULL, NULL, 0}
};

void R_init_RcppClassicExamples(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
