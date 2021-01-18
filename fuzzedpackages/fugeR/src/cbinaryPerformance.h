#ifndef _fugeR_CBINARYPERFORMANCE_H
#define _fugeR_CBINARYPERFORMANCE_H

#include <Rcpp.h>

RcppExport SEXP cbinaryPerformance(  SEXP lstPredicted,
									 SEXP lstActual,
									 SEXP threshold ) ;

#endif
