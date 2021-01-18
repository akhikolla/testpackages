// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdlib.h> // for NULL

// all MixAll method that can be used from R are defined there
#include "../inst/include/MixAll.h"

// declare functions
static const R_CallMethodDef callMethods[]  =
{
  {"clusterMixture",    (DL_FUNC) &clusterMixture, 4},
  {"clusterMixedData",  (DL_FUNC) &clusterMixedData, 3},
  {"clusterPredict",    (DL_FUNC) &clusterPredict, 3},
  {"computeGramMatrix", (DL_FUNC) &computeGramMatrix, 3},
  {"kmm",               (DL_FUNC) &kmm, 4},
  {"kmmMixedData",      (DL_FUNC) &kmmMixedData, 3},
  {"learnMixture",      (DL_FUNC) &learnMixture, 4},
  {"learnMixedData",    (DL_FUNC) &learnMixedData, 3},
  {"learnKmm",          (DL_FUNC) &learnKmm, 4},
//  {"learnKmm",          (DL_FUNC) &learnKmmMixedData, 3},
  {NULL}
};


void R_init_myRoutines(DllInfo *info)
{
	/* Register the .Call routines.
	No .C  .Fortran() or .External() routines,
	so pass those arrays as NULL.
	*/
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
}
