#ifndef _fugeR_CDEFUZZIFY_H
#define _fugeR_CDEFUZZIFY_H

#include <Rcpp.h>

RcppExport SEXP cdefuzzify(SEXP nbCase,
						 SEXP nbRule,
						 SEXP inRule,
						 SEXP lstMf,
						 SEXP lstMfId,
						 SEXP defaultMfId,
						 SEXP lstActivation ) ;

#endif
