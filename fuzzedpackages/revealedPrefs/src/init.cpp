#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// Register routines, to avoid CRAN warning
// "Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'"

RcppExport SEXP CheckWarp(SEXP x, SEXP p, SEXP afriat);
RcppExport SEXP CheckSarp(SEXP x, SEXP p, SEXP afriat);
RcppExport SEXP DeepSarp(SEXP quanti, SEXP prices, SEXP afriat);
RcppExport SEXP CheckGarp(SEXP px, SEXP afriat);
RcppExport SEXP CpUp(SEXP px, SEXP samples, SEXP afriat);
RcppExport SEXP FastUp(SEXP px, SEXP samples, SEXP afriat);
RcppExport SEXP CpLow(SEXP x, SEXP p, SEXP samples, SEXP afriat);
RcppExport SEXP IndirectPrefs(SEXP px, SEXP afriat);
RcppExport SEXP SimAxiom(SEXP nobs, SEXP ngoods, SEXP afriat, SEXP maxit, 
                         SEXP pmin, SEXP pmax, SEXP qmin, SEXP qmax, 
                         SEXP axiom);

static const R_CallMethodDef callMethods[] = {
  {"CheckWarp", (DL_FUNC) &CheckWarp, 3},
  {"CheckSarp", (DL_FUNC) &CheckSarp, 3},
  {"DeepSarp", (DL_FUNC) &DeepSarp, 3},
  {"CheckGarp", (DL_FUNC) &CheckGarp, 2},
  {"CpUp", (DL_FUNC) &CpUp, 3},
  {"FastUp", (DL_FUNC) &FastUp, 3},
  {"CpLow", (DL_FUNC) &CpLow, 4},
  {"IndirectPrefs", (DL_FUNC) &IndirectPrefs, 2},
  {"SimAxiom", (DL_FUNC) &SimAxiom, 9},
  NULL
};

void R_init_revealedPrefs(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
