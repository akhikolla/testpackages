// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
#include "cblas.h"

typedef vec (*func)(vec x, void* arg, int n);
typedef void* (*initfunc)(List arg);
typedef void (*freefunc)(void *arg);

struct weib {
  int nobs; // number of observations
  int natt; // number of attributes
  int nsmp; // number of MCMC samples
  //double *X; // row-major
  mat X; // nobs x natt
  vec alpha; // nsmp
  mat beta; // natt x nsmp
  mat Xbeta; // nobs x nsmp; used in optimized method
};

vec weib_sfunc(vec t, void *arg, int n) {
  weib *argc = (weib*)arg;
  int nsmp = argc->nsmp, nobs = argc->nobs, natt = argc->natt;
  int idx_smp = n / nobs;
  int idx_obs = n - idx_smp * nobs;

  mat X = argc->X; // make sure we are not copying data each time
  mat beta = argc->beta;
  vec alpha = argc->alpha;
  
  double exbeta = exp(accu(X.row(idx_obs) % beta.row(idx_smp)));
  return exp(- pow(t, alpha(idx_smp)) * exbeta);
}
// [[Rcpp::export]]
XPtr<func> weib_getPtr_func() {
  XPtr<func> p(new func(weib_sfunc), true);
  return p;
}

void* weib_init(List arg) { // add an optimized version that pre-computes Xbeta
  weib* myweib = new weib;
  myweib->nobs = arg[0];
  myweib->natt = arg[1];
  myweib->nsmp = arg[2];
  myweib->X = mat(REAL(arg[3]), myweib->nobs, myweib->natt, true, true);
  myweib->alpha = vec(REAL(arg[4]), myweib->nsmp, true, true);
  myweib->beta = mat(REAL(arg[5]), myweib->nsmp, myweib->natt, true, true);
  return (void*)myweib;
}
// [[Rcpp::export]]
XPtr<initfunc> weib_getPtr_init() {
  XPtr<initfunc> p(new initfunc(weib_init), true);
  return p;
}

void weib_free(void *arg) {
  delete (weib*)arg;
}
// [[Rcpp::export]]
XPtr<freefunc> weib_getPtr_free() {
  XPtr<freefunc> p(new freefunc(weib_free), true);
  return p;
}
