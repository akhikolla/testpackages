#ifndef _SpatialTools_GENERATE_NORM_H
#define _SpatialTools_GENERATE_NORM_H

#include <RcppArmadillo.h>

RcppExport SEXP rmvnorm(SEXP nsims, SEXP mus, SEXP Vs, SEXP methods);

RcppExport SEXP condnorm_par(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP coeffs, SEXP Xs, SEXP Xps, SEXP methods);

#endif
