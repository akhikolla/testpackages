#ifndef _bayesImageS_SMC_POTTS_H
#define _bayesImageS_SMC_POTTS_H

#include <RcppArmadillo.h>

RcppExport SEXP smcPotts(SEXP y, SEXP neighbors, SEXP blocks, SEXP param, SEXP priors);

RcppExport SEXP initSedki(SEXP y, SEXP neighbors, SEXP blocks, SEXP param, SEXP priors);

RcppExport SEXP testResample(SEXP values, SEXP weights, SEXP pseudo);

RcppExport SEXP exactPotts(SEXP neighbors, SEXP blocks, SEXP k, SEXP beta);

#endif
