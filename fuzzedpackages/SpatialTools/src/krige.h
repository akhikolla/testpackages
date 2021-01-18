#ifndef _SpatialTools_KRIGE_H
#define _SpatialTools_KRIGE_H

#include <RcppArmadillo.h>

RcppExport SEXP pweights_uk(SEXP Xs, SEXP Vs, SEXP Xps, SEXP Vps, SEXP Vops);

RcppExport SEXP mspe_uk(SEXP ws, SEXP Vs, SEXP Vps, SEXP Vops);

RcppExport SEXP krige_uk(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP Xs, SEXP Xps, SEXP nsims, 
	SEXP Vediags, SEXP methods);

RcppExport SEXP krige_ok(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP nsims, SEXP Vediags, SEXP methods);

RcppExport SEXP krige_sk(SEXP ys, SEXP Vs, SEXP Vps, SEXP Vops, SEXP ms, SEXP nsims, SEXP Vediags, SEXP methods);

RcppExport SEXP spLMPredict(SEXP ys, SEXP coordss, SEXP pcoordss,
                            SEXP Xs, SEXP Xps,
                            SEXP Bs,
                            SEXP sigmasqs, SEXP phis, SEXP nus,
                            SEXP evs, SEXP fvs,
                            SEXP cov_models,
                            SEXP methods, SEXP nreports, SEXP verboses);
#endif
