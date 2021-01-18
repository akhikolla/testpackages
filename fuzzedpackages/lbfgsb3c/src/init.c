#include "../inst/include/lbfgsb3c.h"
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

void lbfgsb3C_(int n, int lmm, double *x, double *lower,
	       double *upper, int *nbd, double *Fmin, optimfn fn,
	       optimgr gr, int *fail, void *ex, double factr,
	       double pgtol, int *fncount, int *grcount,
	       int maxit, char *msg, int trace, int iprint,
	       double atol, double rtol, double *g);

SEXP _lbfgsb3c_lbfgsb3cpp(SEXP, SEXP, SEXP, SEXP, SEXP,
			  SEXP, SEXP);

void R_init_lbfgsb3c(DllInfo *info){
  R_CallMethodDef callMethods[]  = {
    {"_lbfgsb3c_lbfgsb3cpp", (DL_FUNC) &_lbfgsb3c_lbfgsb3cpp, 7},
    {NULL, NULL, 0}
  };
  // C callable to assign environments.
  R_RegisterCCallable("lbfgsb3c", "lbfgsb3C_", (DL_FUNC) lbfgsb3C_);
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
