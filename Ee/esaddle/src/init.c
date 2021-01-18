#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP ecgfCpp(SEXP lambda_, SEXP X_, SEXP mix_, SEXP grad_, SEXP kum1_, SEXP kum2_);

static const R_CallMethodDef CallEntries[] = {
  {"ecgfCpp", (DL_FUNC) &ecgfCpp, 6},
  {NULL, NULL, 0}
};

void R_init_esaddle(DllInfo *info)
{
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
