#include <R_ext/Rdynload.h>
#include "mcmcPotts.h"
#include "smcPotts.h"
#include "PottsUtil.h"

static const R_CallMethodDef callEntries[]  = {
  {"mcmcPotts", (DL_FUNC) &mcmcPotts, 8},
  {"mcmcPottsNoData", (DL_FUNC) &mcmcPottsNoData, 6},
  {"swNoData", (DL_FUNC) &swNoData, 6},
  {"gibbsGMM", (DL_FUNC) &gibbsGMM, 4},
  {"gibbsNorm", (DL_FUNC) &gibbsNorm, 3},
  {"sufficientStat", (DL_FUNC) &sufficientStat, 4},
  {"smcPotts", (DL_FUNC) &smcPotts, 5},
  {"initSedki", (DL_FUNC) &initSedki, 5},
  {"testResample", (DL_FUNC) &testResample, 3},
  {"exactPotts", (DL_FUNC) &exactPotts, 4},
  NULL
};

void R_init_bayesImageS(DllInfo *info) {
  R_registerRoutines(info,
                     NULL, callEntries,
                     NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
