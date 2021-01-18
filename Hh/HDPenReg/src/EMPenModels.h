#ifndef EMPENMODELS_H
#define EMPENMODELS_H

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

SEXP EMlassoMain( SEXP data, SEXP response
                 , SEXP lambda, SEXP intercept
                 , SEXP maxStep, SEXP burn
                 , SEXP threshold, SEXP eps, SEXP epsCG);

SEXP EMlogisticLassoMain( SEXP data, SEXP response
                         , SEXP lambda, SEXP intercept
                         , SEXP maxStep, SEXP burn
                         , SEXP threshold, SEXP eps, SEXP epsCG);

SEXP EMfusedLassoMain( SEXP data, SEXP response
                      , SEXP lambda1, SEXP lambda2, SEXP intercept
                      , SEXP maxStep, SEXP burn
                      , SEXP eps, SEXP eps0, SEXP epsCG);

SEXP EMlogisticFusedLassoMain( SEXP data, SEXP response
                              , SEXP lambda1, SEXP lambda2, SEXP intercept
                              , SEXP maxStep, SEXP burn
                              , SEXP eps, SEXP eps0, SEXP epsCG);

SEXP cvEMlassoMain( SEXP data, SEXP response
                   , SEXP lambda
                   , SEXP nbFolds
                   , SEXP intercept
                   , SEXP maxStep, SEXP burn
                   , SEXP threshold, SEXP eps, SEXP epsCG);

SEXP cvEMfusedLasso1DMain( SEXP data, SEXP response
                          , SEXP lambda1, SEXP lambda2
                          , SEXP optimL1, SEXP nbFolds
                          , SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

SEXP cvEMfusedLasso2DMain( SEXP data, SEXP response
                          , SEXP lambda1, SEXP lambda2
                          , SEXP nbFolds
                          , SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

SEXP cvEMlogisticLassoMain( SEXP data, SEXP response
                           , SEXP lambda, SEXP nbFolds
                           , SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

SEXP cvEMlogisticFusedLasso1DMain( SEXP data, SEXP response
                                  , SEXP lambda1, SEXP lambda2
                                  , SEXP optimL1, SEXP nbFolds
                                  , SEXP intercept, SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

SEXP cvEMlogisticFusedLasso2DMain( SEXP data, SEXP response
                                  , SEXP lambda1, SEXP lambda2
                                  , SEXP nbFolds
                                  , SEXP intercept
                                  , SEXP maxStep, SEXP burn, SEXP threshold, SEXP eps, SEXP epsCG);

#ifdef __cplusplus
}
#endif /*__cplusplus */


#endif
