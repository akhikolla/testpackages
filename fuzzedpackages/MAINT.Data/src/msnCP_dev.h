#ifndef _msnCP_dev_H
#define _msnCP_dev_H

#include "RcppArmadillo.h"

using namespace Rcpp ;
using namespace arma ;

const double b = sqrt(2./PI);

void Rprintivctzzz(const int p,const IntegerVector& v);

template<class VT>
void Rprintv(const int p,const VT& v)
{
  for (int i=0;i<p;i++)  Rprintf("%f ",v(i));
  Rprintf("\n");
} 

template<class MT>
void RprintM(const int m,const int n,const MT& M)
{
  for (int r=0;r<m;r++) {
    for (int c=0;c<n;c++) Rprintf("%f ",M(r,c));
    Rprintf("\n");
  }
  Rprintf("\n");
}

template<class VT>
void Rprintv_rctb(const int p,const VT& v)
{
  for (int i=0;i<p;i++)  Rprintf("%f ",v[i]);
  Rprintf("\n");
} 

template<class MT>
void RprintM_rctb(const int m,const int n,const MT& M)
{
  for (int r=0;r<m;r++) {
    for (int c=0;c<n;c++) Rprintf("%f ",M[r,c]);
    Rprintf("\n");
  }
  Rprintf("\n");
}

inline int ncovp(int Config,int q,int p)
{
  switch (Config)  {
    case 1: return p*(p+1)/2;
    case 2: return p+q+q*(q-1);
    case 3: return p+q;
    case 4: return p+q*(q-1);
    case 5: return p;
  }
  return 0;
}

void cov2cor(const int p,const mat& S,mat& R);
void cnvCPtoDP(const int p,const NumericVector mu,const mat& Sigma,const NumericVector gamma1,
		vec& ksi, mat& Omega, vec& alpha, mat& Omegabar, vec& delta,
		double* c2, bool* admissible, double* viol, const double c2tol, const double ldRtol, const double limlnk2, const bool FixedArrays);
double msnCP_dev1(NumericVector& param, const NumericMatrix& y, const IntegerVector& grpind, 
		const int Config, const int n, const int p, const int k, const double limlnk2, 
		const bool trace, const double c2tol, const double ldRtol, 
		const double PenF, const double PenC, const bool nopenalty,
		const double MachineEPS, const bool FixedArrays, const bool Srpar);


RcppExport SEXP msnCP_dev(SEXP param_s, SEXP y_s, SEXP grpind_s, SEXP Config_s, SEXP n_s, SEXP p_s, SEXP k_s, SEXP limlnk2_s, 
  SEXP trace_s, SEXP c2tol_s, SEXP ldRtol_s, SEXP PenF_s, SEXP PenC_s, SEXP nopenalty_s, SEXP MachineEPS_s, SEXP Srpar_s);

#endif
