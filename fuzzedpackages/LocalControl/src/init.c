#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _LocalControl_getMaxDist(SEXP);
extern SEXP _LocalControl_newCRLC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _LocalControl_newLC(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_LocalControl_getMaxDist", (DL_FUNC) &_LocalControl_getMaxDist, 1},
    {"_LocalControl_newCRLC",    (DL_FUNC) &_LocalControl_newCRLC,    5},
    {"_LocalControl_newLC",      (DL_FUNC) &_LocalControl_newLC,      3},
    {NULL, NULL, 0}
};

void R_init_LocalControl(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

