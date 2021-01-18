#ifndef _LARSR_H
#define _LARSR_H

#ifdef _OPENMP
#include <omp.h>
#endif


#ifdef __cplusplus
extern "C"
{
#endif

SEXP larsmain(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps);
SEXP cvlarsmain(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps, SEXP nbFold, SEXP partition, SEXP index, SEXP mode);
SEXP fusionmain(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps);

#ifdef __cplusplus
}
#endif /*__cplusplus */


#endif
