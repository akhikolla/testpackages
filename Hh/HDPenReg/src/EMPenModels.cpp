#include <Rcpp.h>
#include "EMPenModels.h"



using namespace Rcpp;

// lasso and logistic lasso
RcppExport SEXP EMlassoC( SEXP data, SEXP response
                       , SEXP lambda, SEXP intercept
                       , SEXP maxStep, SEXP burn
                       , SEXP threshold, SEXP eps, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(1);
#endif
  return EMlassoMain(data, response, lambda, intercept, maxStep, burn, threshold, eps, epsCG);
}
RcppExport SEXP EMlogisticLasso( SEXP data, SEXP response
                               , SEXP lambda, SEXP intercept
                               , SEXP maxStep, SEXP burn
                               , SEXP threshold, SEXP eps, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(1);
#endif
  return EMlogisticLassoMain(data, response, lambda, intercept, maxStep, burn, threshold, eps, epsCG);
}

// fused lasso and logistic fused lasso
RcppExport SEXP EMfusedLasso( SEXP data, SEXP response
                            , SEXP lambda1, SEXP lambda2, SEXP intercept
                            , SEXP maxStep, SEXP burn
                            , SEXP eps, SEXP eps0, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(1);
#endif
  return EMfusedLassoMain(data, response, lambda1, lambda2, intercept, maxStep, burn, eps, eps0, epsCG);
}

RcppExport SEXP EMlogisticFusedLasso( SEXP data, SEXP response
                                    , SEXP lambda1, SEXP lambda2, SEXP intercept
                                    , SEXP maxStep, SEXP burn
                                    , SEXP eps, SEXP eps0, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(1);
#endif
  return EMlogisticFusedLassoMain(data, response, lambda1, lambda2, intercept, maxStep, burn, eps, eps0, epsCG);
}

// cv models
RcppExport SEXP cvEMlasso( SEXP data, SEXP response
                         , SEXP lambda, SEXP nbFolds, SEXP intercept
                         , SEXP maxStep, SEXP burn
                         , SEXP threshold, SEXP eps, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(2);
#endif
  return cvEMlassoMain(data, response, lambda, nbFolds, intercept, maxStep, burn, threshold, eps, epsCG);
}

RcppExport SEXP cvEMfusedLasso1D( SEXP data, SEXP response
                               , SEXP lambda1, SEXP lambda2
                               , SEXP optimL1, SEXP nbFolds
                               , SEXP intercept, SEXP maxStep
                               , SEXP burn, SEXP threshold
                               , SEXP eps, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(2);
#endif
  return cvEMfusedLasso1DMain( data,  response,  lambda1,  lambda2,  optimL1,  nbFolds,  intercept,  maxStep,  burn,  threshold,  eps,  epsCG);
}

RcppExport SEXP cvEMfusedLasso2D(SEXP data, SEXP response
                               , SEXP lambda1, SEXP lambda2
                               , SEXP nbFolds, SEXP intercept
                               , SEXP maxStep, SEXP burn
                               , SEXP threshold, SEXP eps, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(2);
#endif
  return cvEMfusedLasso2DMain( data,  response, lambda1,  lambda2,  nbFolds,  intercept,  maxStep,  burn,  threshold,  eps,  epsCG);
}

// duplication pour cv logistic, trouver autre chose de mieux
//
RcppExport SEXP cvEMlogisticLasso( SEXP data, SEXP response
                                  , SEXP lambda, SEXP nbFolds
                                  , SEXP intercept, SEXP maxStep
                                  , SEXP burn, SEXP threshold
                                  , SEXP eps, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(2);
#endif
  return cvEMlogisticLassoMain(data, response, lambda, nbFolds, intercept, maxStep, burn, threshold, eps, epsCG);
}

//
RcppExport SEXP cvEMlogisticFusedLasso1D( SEXP data, SEXP response
                                         , SEXP lambda1, SEXP lambda2
                                         , SEXP optimL1, SEXP nbFolds
                                         , SEXP intercept, SEXP maxStep
                                         , SEXP burn, SEXP threshold
                                         , SEXP eps, SEXP epsCG)
{
#ifdef _OPENMP
omp_set_num_threads(2);
#endif
  return cvEMlogisticFusedLasso1DMain( data,  response,  lambda1,  lambda2,  optimL1,  nbFolds,  intercept,  maxStep,  burn,  threshold,  eps,  epsCG);
}

//
RcppExport SEXP cvEMlogisticFusedLasso2D( SEXP data, SEXP response
                                            , SEXP lambda1, SEXP lambda2
                                            , SEXP nbFolds
                                            , SEXP intercept
                                            , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG)
{
#ifdef _OPENMP
  omp_set_num_threads(2);
#endif
  return cvEMlogisticFusedLasso2DMain(data,  response,  lambda1,  lambda2,  nbFolds,  intercept,  maxStep,  burn,  threshold,  eps,  epsCG);
}
