#include <RcppArmadillo.h>

RcppExport SEXP CCK(SEXP Y, SEXP k);
RcppExport SEXP PVCk(SEXP Y, SEXP k);
RcppExport SEXP TIK(SEXP Y, SEXP U, SEXP k, SEXP method);
RcppExport SEXP TIKlc(SEXP Y, SEXP U, SEXP k, SEXP method);
RcppExport SEXP TIK1(SEXP Y, SEXP U, SEXP k);
RcppExport SEXP TSAVE(SEXP X, SEXP slices, SEXP h, SEXP k);
RcppExport SEXP TSIR(SEXP X, SEXP slices, SEXP h, SEXP k);
RcppExport SEXP lblinM(SEXP X, SEXP k);
RcppExport SEXP lbsqM(SEXP X, SEXP k);
RcppExport SEXP EIGEN(SEXP X);
RcppExport SEXP PREPBSS(SEXP X, SEXP n);
