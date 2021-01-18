#ifndef _SpatialTools_COV_TOOLS_H
#define _SpatialTools_COV_TOOLS_H

#include <RcppArmadillo.h>

RcppExport SEXP decomp_cov(SEXP Vs, SEXP methods);

//RcppExport SEXP kweights_uk_arma(SEXP Xs, SEXP Vs, SEXP Xps, SEXP Vps, SEXP Vops);

//RcppExport SEXP mspe_uk_arma(SEXP ws, SEXP Vs, SEXP Vps, SEXP Vops);

#endif
