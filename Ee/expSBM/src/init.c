#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _expSBM_cpp_expSBM_ELBO(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _expSBM_cpp_expSBM_EM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_expSBM_cpp_expSBM_ELBO", (DL_FUNC) &_expSBM_cpp_expSBM_ELBO,  9},
  {"_expSBM_cpp_expSBM_EM",   (DL_FUNC) &_expSBM_cpp_expSBM_EM,   11},
  {NULL, NULL, 0}
};

void R_init_expSBM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


/* TO GENERATE THIS FILE I HAVE USED:
   tools::package_native_routine_registration_skeleton(".")

   OR
   tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
*/
