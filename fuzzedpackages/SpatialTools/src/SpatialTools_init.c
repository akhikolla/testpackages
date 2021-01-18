#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .C calls */
  extern void dist1_c(void *, void *, void *, void *);
extern void dist2_c(void *, void *, void *, void *, void *, void *);

/* .Call calls */
  extern SEXP coincident_cpp(SEXP, SEXP, SEXP);
extern SEXP krige_ok(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP krige_sk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP spLMPredict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"dist1_c", (DL_FUNC) &dist1_c, 4},
  {"dist2_c", (DL_FUNC) &dist2_c, 6},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"coincident_cpp", (DL_FUNC) &coincident_cpp,  3},
  {"krige_ok",       (DL_FUNC) &krige_ok,        7},
  {"krige_sk",       (DL_FUNC) &krige_sk,        8},
  {"spLMPredict",    (DL_FUNC) &spLMPredict,    15},
  {NULL, NULL, 0}
};

void R_init_SpatialTools(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
