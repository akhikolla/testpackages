#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include "kergp.h"

#include <R_ext/Rdynload.h>

static R_NativePrimitiveArgType covm_t[] = { 
  REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType covmm_t[] = { 
  REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType kern1_t[] = { 
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP
};

static const R_CMethodDef CEntries[] = {
  {"C_covWhiteNoise",  (DL_FUNC) &C_covWhiteNoise, 8, covm_t}, 
  {"C_covMat1Mat2_WhiteNoise",  (DL_FUNC) &C_covMat1Mat2_WhiteNoise, 7, covmm_t},
  {"kern1Gauss", (DL_FUNC) &kern1Gauss, 6, kern1_t},
  {"kern1Exp", (DL_FUNC) &kern1Exp, 6, kern1_t},
  {"kern1PowExp", (DL_FUNC) &kern1PowExp, 6, kern1_t},
  {NULL, NULL, 0}
};

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
  CALLDEF(scores_covMan, 5),
  CALLDEF(covMat_covMan, 6),
  CALLDEF(covMatMat_covMan, 7), 
  CALLDEF(varVec_covMan, 6), 
  CALLDEF(covMat_covTS, 7),
  CALLDEF(covMatMat_covTS, 8), 
  CALLDEF(varVec_covTS, 7),
  CALLDEF(scores_covTS, 6),
  CALLDEF(k1ExpC, 3),
  CALLDEF(k1GaussC, 3),  
  CALLDEF(k1PowExpC, 3),
  CALLDEF(k1Matern3_2C, 3),
  CALLDEF(k1Matern5_2C, 3),
  CALLDEF(corLev_CompSymm, 4),
  CALLDEF(corLev_Symm, 4),
  CALLDEF(corLev_LowRank, 5),
  CALLDEF(k1FunExpC, 1),
  CALLDEF(k1FunMatern3_2C, 1),
  CALLDEF(k1FunMatern5_2C, 1), 
  CALLDEF(k1FunGaussC, 1),
  CALLDEF(k1FunPowExpC, 2),  
  {NULL, NULL, 0}
};


void R_init_kergp(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
