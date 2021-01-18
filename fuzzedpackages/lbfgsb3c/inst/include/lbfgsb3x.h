#ifndef __LBFGSB3_H__
#define __LBFGSB3_H__
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

typedef double optimfn(int n, double *par, void *ex);

typedef void optimgr(int n, double *par, double *gr, void *ex);

typedef void (*lbfgsb3_fn)(int n, int lmm, double *x, double *lower,
			   double *upper, int *nbd, double *Fmin, optimfn fn,
			   optimgr gr, int *fail, void *ex, double factr,
			   double pgtol, int *fncount, int *grcount,
			   int maxit, char *msg, int trace, int nREPORT, double atol, double rtol,
			   double *g);

void lbfgsb3C(int n, int lmm, double *x, double *lower,
	      double *upper, int *nbd, double *Fmin, optimfn fn,
	      optimgr gr, int *fail, void *ex, double factr,
	      double pgtol, int *fncount, int *grcount,
	      int maxit, char *msg, int trace, int nREPORT, double atol, double rtol,
	      double *g){
  static lbfgsb3_fn fun=NULL;
  if (fun == NULL) fun = (lbfgsb3_fn) R_GetCCallable("lbfgsb3x","lbfgsb3C_");
  fun(n, lmm, x, lower, upper, nbd, Fmin, fn, gr, fail, ex, factr, pgtol, fncount,
      grcount, maxit, msg, trace, nREPORT, atol, rtol, g);
}
#endif
