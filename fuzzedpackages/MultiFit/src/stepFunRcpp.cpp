#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

typedef struct {
  double ylow;
  double yhigh;
  double f1;
  double f2;
  int kind;
} appr_meth;
 
static double approx1(double v, double *x, double *y, int n,
                      appr_meth *Meth)
{
  int i, j, ij;
  if(!n) return R_NaN;
  i = 0;
  j = n - 1;

  /* handle out-of-domain points */
  if(v < x[i]) return Meth->ylow;
  if(v > x[j]) return Meth->yhigh;

  /* find the correct interval by bisection */
  while(i < j - 1) { /* x[i] <= v <= x[j] */
    ij = (i + j)/2; /* i+1 <= ij <= j-1 */
    if(v < x[ij]) j = ij; else i = ij;
    /* still i < j */
  }
  /* interpolation */
  if(v == x[j]) return y[j];
  if(v == x[i]) return y[i];
  if(Meth->kind == 1) /* linear */
  return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
  else /* 2 : constant */
  return (Meth->f1 != 0.0 ? y[i] * Meth->f1 : 0.0)
    + (Meth->f2 != 0.0 ? y[j] * Meth->f2 : 0.0);
}/* approx1() */
   
static void
  R_approxfun(double *x, double *y, int nxy, double *xout, double *yout,
              int nout, int method, double yleft, double yright, double f)
  {
    int i;
    appr_meth M = {0.0, 0.0, 0.0, 0.0, 0}; /* -Wall */

  M.f2 = f;
  M.f1 = 1 - f;
  M.kind = method;
  M.ylow = yleft;
  M.yhigh = yright;
  for(i = 0; i < nout; i++)
    yout[i] = ISNAN(xout[i]) ? xout[i] : approx1(xout[i], x, y, nxy, &M);
  }
 
NumericVector Approx(NumericVector x, NumericVector y, NumericVector v, 
                     int method, double yleft, double yright, double f) {
  NumericVector xout = v;
  int nx = x.length(), nout = xout.length();
  NumericVector yout = v;
  R_approxfun(REAL(x), REAL(y), nx, REAL(xout), REAL(yout), nout,
              method, yleft, yright, f);
  return yout;
}

NumericVector stepFun(NumericVector z, NumericVector v) 
{
    IntegerVector idx = seq_along(z);
    NumericVector y = z[idx-1];
    NumericVector x = y;
    return Approx(x, y, v, 2, 0, y(y.length()-1), 0);
}

//[[Rcpp::export]]
NumericVector genStepFunCpp( Rcpp::List SPS, NumericVector sort_p, int lc, int lp ) {
  
  NumericVector pp(lc);
  pp.fill(0);
  
  for (int s=0; s<lp; s++) {
    NumericVector sps = Rcpp::NumericVector(SPS(s));
    NumericVector foo;
    if (s==0) {
      foo = sort_p[0];
    } else {
      foo = sort_p[seq(0, std::min(s,lc-1))];
    }
    NumericVector tp = stepFun(sps, foo);
    if (s==0) {
      pp(0) = Rcpp::as<double>(tp);
    } else {
      pp[seq(0, tp.length()-1)] += tp;
    }
  }
  return(pp);
}
