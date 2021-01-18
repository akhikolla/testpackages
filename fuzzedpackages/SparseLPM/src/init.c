#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP SparseLPM_cpp_SLPM_ELBO(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SparseLPM_cpp_SLPM_Optimisation(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"SparseLPM_cpp_SLPM_ELBO",         (DL_FUNC) &SparseLPM_cpp_SLPM_ELBO,         13},
    {"SparseLPM_cpp_SLPM_Optimisation", (DL_FUNC) &SparseLPM_cpp_SLPM_Optimisation, 18},
    {NULL, NULL, 0}
};

void R_init_SparseLPM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

/* TO GENERATE THIS FILE I HAVE USED:
   tools::package_native_routine_registration_skeleton(".")

   OR
   tools::package_native_routine_registration_skeleton(".", character_only = FALSE)
*/
