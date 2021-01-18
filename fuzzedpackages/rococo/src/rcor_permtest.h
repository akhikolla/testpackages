#ifndef _RCOR_PERMTEST_H
#define _RCOR_PERMTEST_H

#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP rcor_permtest(SEXP matx, SEXP maty, SEXP tnorm, SEXP tests,
			      SEXP ogamma, SEXP alternative,
			      SEXP storeValues);

#endif
