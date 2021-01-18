#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <cmath>
#include "Rdtq_types.h"

double callViaXPtr(const double x, SEXP xpsexp)
{
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  double y = fun(x);
  return(y);
}

// gaussian pdf, non-vectorized version of GSL function
static inline double gaussian_pdf(const double x, const double mu, const double sigma)
{
  double u = (x - mu) / sigma;
  double p = (1 / (sqrt (2.0*M_PI) * sigma)) * exp (-(u*u) / 2.0);
  return p;
}

// log gaussian pdf, non-vectorized version of GSL function
static inline double log_gaussian_pdf(const double x, const double mu, const double sigma)
{
  double y = (x - mu) / sigma;
  double out = -(y*y) / 2.0;
  out -= 0.5*(log(2.0*M_PI));
  out -= log(sigma);
  return out;
}

// [[Rcpp::export]]
List rdtq(double h, double k, int bigm, NumericVector init, double T, SEXP drift, SEXP diffusion, double thresh)
{
  XPtr<funcPtr> driftF(drift);
  XPtr<funcPtr> diffF(diffusion);
  funcPtr driftfun = *driftF;
  funcPtr difffun = *diffF;

  unsigned int veclen = 2*bigm+1;
  int iveclen = 2*bigm+1;

  double h12 = sqrt(h);

  NumericVector oldphatn(veclen);
  NumericVector phatn(veclen);
  NumericVector phatnp1(veclen);
  NumericVector xvec(veclen);
  for (int i=-bigm;i<=bigm;i++) {
    xvec(i+bigm) = i*k;
  }

  int startn;

  if (init.size()==1)
  {
    double initval = init(0);
    // pdf after one time step
    // need to do this if initial condition is a fixed constant
    // in which case initial PDF is a Dirac delta, so we do one step manually
    double mymu = initval + driftfun(initval)*h;
    double mysig = (std::fabs(difffun(initval)))*h12;
    for (int i=0;i<iveclen;i++) {
      phatn(i) = gaussian_pdf(xvec(i),mymu,mysig);
    }
    startn = 1;
  }
  else
  {
    for (int i=0;i<iveclen;i++) phatn(i) = init(i);
    startn = 0;
  }

  NumericVector mymuvec(veclen);
  NumericVector mysigvec(veclen);

  for (int j=-bigm;j<=bigm;j++)
  {
    mymuvec(j+bigm) = j*k + driftfun(j*k)*h;
    mysigvec(j+bigm) = (std::fabs(difffun(j*k)))*h12;
  }

  // iterate
  int bign = ceil(T/h);
  for (int n=startn;n<bign;n++) {
    for (int i=-bigm;i<=bigm;i++) {
      double explogterm = 0.0;
      //double term = 0.0;
      for (int j=-bigm;j<=bigm;j++) {
        //term += gaussian_pdf(i*k,mymuvec(j+bigm),mysigvec(j+bigm)) * phatn(j+bigm);
        if (phatn(j+bigm) < thresh) continue;
        double logterm = log_gaussian_pdf(i*k,mymuvec(j+bigm),mysigvec(j+bigm));
        logterm += log(phatn(j+bigm));
        explogterm += exp(logterm);
      }
      phatnp1(i+bigm) = k*explogterm;
    }
    for (int i=0; i<veclen; i++) phatn(i) = phatnp1(i);
  }

  return List::create(Named("xvec")=xvec,Named("pdf")=phatn);
}

// version of rdtq in which user specifies the boundaries [a,b] of the spatial grid,
// together with the total number of grid points

// [[Rcpp::export]]
List rdtqgrid(double h, double a, double b, unsigned int veclen, NumericVector init, double T, SEXP drift, SEXP diffusion, double thresh)
{
  XPtr<funcPtr> driftF(drift);
  XPtr<funcPtr> diffF(diffusion);
  funcPtr driftfun = *driftF;
  funcPtr difffun = *diffF;

  int iveclen = (int) veclen;

  double h12 = sqrt(h);

  NumericVector oldphatn(veclen);
  NumericVector phatn(veclen);
  NumericVector phatnp1(veclen);
  NumericVector xvec(veclen);
  double k = (b-a)/(iveclen-1);
  for (int i=0;i<(iveclen-1);i++) {
    xvec(i) = a + i*k;
  }
  xvec(veclen-1) = b;

  int startn;

  if (init.size()==1)
  {
    double initval = init(0);
    // pdf after one time step
    // need to do this if initial condition is a fixed constant
    // in which case initial PDF is a Dirac delta, so we do one step manually
    double mymu = initval + driftfun(initval)*h;
    double mysig = (std::fabs(difffun(initval)))*h12;
    for (int i=0;i<iveclen;i++) {
      phatn(i) = gaussian_pdf(xvec(i),mymu,mysig);
    }
    startn = 1;
  }
  else
  {
    for (int i=0;i<iveclen;i++) phatn(i) = init(i);
    startn = 0;
  }

  NumericVector mymuvec(veclen);
  NumericVector mysigvec(veclen);

  for (int j=0;j<veclen;j++)
  {
    mymuvec(j) = xvec(j) + driftfun(xvec(j))*h;
    mysigvec(j) = (std::fabs(difffun(xvec(j))))*h12;
  }

  // iterate
  int bign = ceil(T/h);
  for (int n=startn;n<bign;n++) {
    for (int i=0;i<veclen;i++) {
      double explogterm = 0.0;
      for (int j=0;j<veclen;j++) {
        if (phatn(j) < thresh) continue;
        double logterm = log_gaussian_pdf(xvec(i),mymuvec(j),mysigvec(j));
        logterm += log(phatn(j));
        explogterm += exp(logterm);
      }
      phatnp1(i) = k*explogterm;
    }
    for (int i=0; i<veclen; i++) phatn(i) = phatnp1(i);
  }

  return List::create(Named("xvec")=xvec,Named("pdf")=phatn);
}

