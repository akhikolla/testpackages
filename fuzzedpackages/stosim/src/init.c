/*
package_native_routine_registration_skeleton("C:/Repositories/stosim/pkg/stosim")
*/
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP DetailOpLinesCPP(SEXP, SEXP, SEXP, SEXP);
extern SEXP MultiTrainWithInventoryCPP(SEXP, SEXP, SEXP, SEXP);
extern SEXP SimulationHistory(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"DetailOpLinesCPP",        (DL_FUNC) &DetailOpLinesCPP,         4},
    {"MultiTrainWithInventoryCPP", (DL_FUNC) &MultiTrainWithInventoryCPP,  4},
    {"SimulationHistory",       (DL_FUNC) &SimulationHistory,       12},
    {NULL, NULL, 0}
};

void R_init_stosim(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

