#ifndef _AdaptiveSparsity_RCPP_HELLO_WORLD_H
#define _AdaptiveSparsity_RCPP_HELLO_WORLD_H

#include <RcppArmadillo.h>

RcppExport SEXP olsInit(SEXP inputMatX);
RcppExport SEXP CSL(SEXP it, SEXP init, SEXP k, SEXP n, SEXP eps, SEXP covMat, SEXP ansQ);

#endif
