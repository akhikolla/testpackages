


#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP boundedBinom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP boundedGauss(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP boundedGaussVar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP boundedPoisson(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP confBand(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP forwardBinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP forwardGauss(SEXP, SEXP, SEXP, SEXP);
extern SEXP forwardGaussCut(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP forwardGaussInhibit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP forwardGaussVar(SEXP, SEXP, SEXP);
extern SEXP forwardPoisson(SEXP, SEXP, SEXP);
extern SEXP pathBinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP pathGauss(SEXP, SEXP, SEXP, SEXP);
extern SEXP pathGaussCut(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pathGaussVar(SEXP, SEXP, SEXP);
extern SEXP pathPoisson(SEXP, SEXP, SEXP);
extern SEXP _stepR_callRoutines(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _stepR_colMax(SEXP);
extern SEXP _stepR_criticalValuesWeights(SEXP, SEXP, SEXP);
extern SEXP _stepR_inOrdered(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"boundedBinom",                (DL_FUNC) &boundedBinom,                7},
  {"boundedGauss",                (DL_FUNC) &boundedGauss,                7},
  {"boundedGaussVar",             (DL_FUNC) &boundedGaussVar,             6},
  {"boundedPoisson",              (DL_FUNC) &boundedPoisson,              6},
  {"confBand",                    (DL_FUNC) &confBand,                    6},
  {"forwardBinom",                (DL_FUNC) &forwardBinom,                4},
  {"forwardGauss",                (DL_FUNC) &forwardGauss,                4},
  {"forwardGaussCut",             (DL_FUNC) &forwardGaussCut,             9},
  {"forwardGaussInhibit",         (DL_FUNC) &forwardGaussInhibit,         7},
  {"forwardGaussVar",             (DL_FUNC) &forwardGaussVar,             3},
  {"forwardPoisson",              (DL_FUNC) &forwardPoisson,              3},
  {"pathBinom",                   (DL_FUNC) &pathBinom,                   4},
  {"pathGauss",                   (DL_FUNC) &pathGauss,                   4},
  {"pathGaussCut",                (DL_FUNC) &pathGaussCut,                9},
  {"pathGaussVar",                (DL_FUNC) &pathGaussVar,                3},
  {"pathPoisson",                 (DL_FUNC) &pathPoisson,                 3},
  {"_stepR_callRoutines",          (DL_FUNC) &_stepR_callRoutines,          7},
  {"_stepR_colMax",                (DL_FUNC) &_stepR_colMax,                1},
  {"_stepR_criticalValuesWeights", (DL_FUNC) &_stepR_criticalValuesWeights, 3},
  {"_stepR_inOrdered",             (DL_FUNC) &_stepR_inOrdered,             2},
  {NULL, NULL, 0}
};

void R_init_stepR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
