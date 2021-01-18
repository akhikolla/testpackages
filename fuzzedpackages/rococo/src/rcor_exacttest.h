#ifndef _RCOR_EXACTTEST_H
#define _RCOR_EXACTTEST_H

#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP permNextWrapper(SEXP perm, SEXP sign);

RcppExport SEXP rcor_exacttest(SEXP matx, SEXP maty, SEXP tnorm, SEXP tests,
			       SEXP ogamma, SEXP alternative,
			       SEXP storeValues);

#endif
