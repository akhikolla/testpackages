#include <RcppArmadillo.h>

RcppExport SEXP rnormCube(SEXP varp1, SEXP varp2, SEXP varp3);
RcppExport SEXP eigenVectors(SEXP varx);
RcppExport SEXP symmetricPower(SEXP varx, SEXP varr);
RcppExport SEXP mFOBIMatrix(SEXP varx);
RcppExport SEXP mFOBIMatrixNorm(SEXP varx);
RcppExport SEXP mJADEMatrix(SEXP varx, SEXP vari, SEXP varj, SEXP varcov);
RcppExport SEXP matrixCovariance(SEXP varx);
RcppExport SEXP mAutoCovMatrix(SEXP varx, SEXP varlag);
RcppExport SEXP mTGFOBIMatrix(SEXP varx, SEXP varlag);
RcppExport SEXP mTGJADEMatrix(SEXP varx, SEXP vari, SEXP varj, SEXP varlags);
