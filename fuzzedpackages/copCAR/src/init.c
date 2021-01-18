#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP copCAR_buildM_(SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(pbivnorm)(void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"copCAR_buildM_", (DL_FUNC) &copCAR_buildM_, 3},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"pbivnorm", (DL_FUNC) &F77_NAME(pbivnorm), 7},
    {NULL, NULL, 0}
};

void R_init_copCAR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
