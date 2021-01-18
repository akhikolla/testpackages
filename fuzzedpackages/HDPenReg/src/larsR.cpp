#include <Rcpp.h>
#include "larsR.h"

extern "C" SEXP lars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps)
{ return larsmain(data, response, nbIndiv, nbVar, maxStep, intercept, eps);}

extern "C"  SEXP fusion(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps)
{ return fusionmain(data, response, nbIndiv, nbVar, maxStep, intercept, eps);}

extern "C"  SEXP cvlars(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps, SEXP nbFold, SEXP partition, SEXP index, SEXP mode)
{ return cvlarsmain(data, response, nbIndiv, nbVar, maxStep, intercept, eps, nbFold, partition, index, mode);}
