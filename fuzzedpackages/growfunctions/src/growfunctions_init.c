#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP GP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GPBFIX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GPDPMIX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GPDPMIXMIS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP GPFIX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP IGMRFDPMIX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP IGMRFDPMIXCOUNT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP predict_bb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP predict_gmrf_bb(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"GP",              (DL_FUNC) &GP,              20},
  {"GPBFIX",          (DL_FUNC) &GPBFIX,          21},
  {"GPDPMIX",         (DL_FUNC) &GPDPMIX,         24},
  {"GPDPMIXMIS",      (DL_FUNC) &GPDPMIXMIS,      25},
  {"GPFIX",           (DL_FUNC) &GPFIX,           21},
  {"IGMRFDPMIX",      (DL_FUNC) &IGMRFDPMIX,      21},
  {"IGMRFDPMIXCOUNT", (DL_FUNC) &IGMRFDPMIXCOUNT, 24},
  {"predict_bb",      (DL_FUNC) &predict_bb,       6},
  {"predict_gmrf_bb", (DL_FUNC) &predict_gmrf_bb,  3},
  {NULL, NULL, 0}
};

void R_init_growfunctions(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
