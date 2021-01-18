#ifndef _msnCP_dev_gra_h
#define _msnCP_dev_gra_h

#include "RcppArmadillo.h"
using namespace arma;

void msnCP_ll_grad(const NumericVector& mu1, const NumericMatrix& beta2k, const mat& Sigma, 
  const NumericVector& gamma1, const NumericMatrix& y, const IntegerVector& grpind, 
  const int n, const int p, const int k, const int nvcovpar, const double limlnk2, const bool trace,
  const double PenF, const double ldRtol, const double c2tol, const double beta0tol,
  vec& CPgrad, const double Machinetol, const bool FixedArrays);

NumericVector msnCP_dev_grad1(NumericVector& param, const NumericMatrix& y, const IntegerVector& grpind,
  const int Config, const int n, const int p, const int k, const double limlnk2,
  const bool trace, const double c2tol, const double ldRtol, const double beta0tol, 
  const double PenF, const double MachineEPS, const bool FixedArrays, const bool Srpar);

RcppExport SEXP msnCP_dev_grad(SEXP param_s, SEXP y_s, SEXP grpind_s, SEXP Config_s, SEXP n_s, SEXP p_s, SEXP k_s, SEXP limlnk2_s, 
  SEXP trace_s, SEXP c2tol_s, SEXP ldRtol_s, SEXP beta0tol_s, SEXP PenF_s, SEXP MachineEPS_s, SEXP Srpar_s) ;


#endif
