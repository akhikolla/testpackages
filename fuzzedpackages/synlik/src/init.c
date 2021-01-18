#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP cleanStats(SEXP inMat);
SEXP checkBoundsCpp(SEXP theMean_, SEXP cholFact_, SEXP toCheck_, SEXP upper_, SEXP lower_, SEXP output_);
SEXP simpleModelsWrap(SEXP model, SEXP days, SEXP nSimul, SEXP params, SEXP nBurn, SEXP randInit, SEXP initVal);

void slacf(double *acf,double *x,int *n,int *n_reps,int *max_lag,double *NAcode,int *correlation);
void slnlar(double *beta, double *x,int *n,int *n_reps,int *n_terms,
            int *lag,int *power,double *NAcode);
void order_reg(double *beta, double *x,double *z,int *n,int *n_reps,int *np,int *diff);
void blowC(double *n,double *theta,double *e,double *e1,int *burn_in,int *n_t, int *n_reps);

R_CallMethodDef CallEntries[] = {
  {"cleanStats", (DL_FUNC) &cleanStats, 1},
  {"checkBoundsCpp", (DL_FUNC) &checkBoundsCpp, 6},
  {"simpleModelsWrap", (DL_FUNC) &simpleModelsWrap, 7},
  {NULL, NULL, 0}
};
	
R_CMethodDef cMethods[] = {
   {"slacf", (DL_FUNC) &slacf, 7},
   {"slnlar", (DL_FUNC) &slnlar, 8},
   {"order_reg", (DL_FUNC) &order_reg, 7},
   {"blowC", (DL_FUNC) &blowC, 7},
   {NULL, NULL, 0}
};

void R_init_synlik(DllInfo *info)
{
  R_registerRoutines(info, cMethods, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}


