// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "rtkore.h"

// declare functions
static const R_CallMethodDef callMethods[]  =
{
  {"stk_version", (DL_FUNC) &stk_version, 1},
  {"fastBetaRand", (DL_FUNC) &fastBetaRand, 3},
  {"fastBinomialRand", (DL_FUNC) &fastBinomialRand, 3},
  {"fastCategoricalRand", (DL_FUNC) &fastCategoricalRand,  2},
  {"fastCauchyRand", (DL_FUNC) &fastCauchyRand, 3},
  {"fastChiSquaredRand", (DL_FUNC) &fastChiSquaredRand, 2},
  {"fastExponentialRand", (DL_FUNC) &fastExponentialRand, 2},
  {"fastFisherSnedecorRand", (DL_FUNC) &fastFisherSnedecorRand, 3},
  {"fastGammaRand", (DL_FUNC) &fastGammaRand, 3},
  {"fastGeometricRand", (DL_FUNC) &fastGeometricRand, 2},
  {"fastHyperGeometricRand", (DL_FUNC) &fastHyperGeometricRand, 4},
  {"fastLogisticRand", (DL_FUNC) &fastLogisticRand, 3},
  {"fastLogNormalRand", (DL_FUNC) &fastLogNormalRand, 3},
  {"fastNegativeBinomialRand", (DL_FUNC) &fastNegativeBinomialRand, 3},
  {"fastNormalRand", (DL_FUNC) &fastNormalRand, 3},
  {"fastPoissonRand", (DL_FUNC) &fastPoissonRand, 2},
  {"fastStudentRand", (DL_FUNC) &fastStudentRand, 2},
  {"fastUniformRand", (DL_FUNC) &fastUniformRand, 3},
  {"fastUniformDiscreteRand", (DL_FUNC) &fastUniformDiscreteRand, 3},
  {"fastWeibullRand", (DL_FUNC) &fastWeibullRand, 3},
  {NULL}
};


void R_init_myRoutines(DllInfo *info)
{
	/* Register the .Call routines.
	No .C  .Fortran() or .External() routines,
	so pass those arrays as NULL.
	*/
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}
