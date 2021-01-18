#ifndef _bayesImageS_MCMC_POTTS_H
#define _bayesImageS_MCMC_POTTS_H

#include <RcppArmadillo.h>

RcppExport SEXP mcmcPotts(SEXP y, SEXP neighbors, SEXP blocks, SEXP niter, SEXP nburn, SEXP priors, SEXP mh, SEXP truth);

RcppExport SEXP mcmcPottsNoData(SEXP beta, SEXP k, SEXP neighbors, SEXP blocks, SEXP niter, SEXP random);

RcppExport SEXP swNoData(SEXP beta, SEXP k, SEXP neighbors, SEXP blocks, SEXP niter, SEXP random);

RcppExport SEXP gibbsGMM(SEXP y, SEXP niter, SEXP nburn, SEXP priors);

RcppExport SEXP gibbsNorm(SEXP y, SEXP niter, SEXP priors);

RcppExport SEXP relabel(SEXP y, SEXP neighbors, SEXP blocks, SEXP samples, SEXP truth);

RcppExport SEXP sufficientStat(SEXP zMx, SEXP neighbors, SEXP blocks, SEXP k);

#endif
