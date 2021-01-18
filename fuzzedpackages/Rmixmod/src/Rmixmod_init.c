#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP clusteringMain(SEXP);
extern SEXP learnMain(SEXP);
extern SEXP predictMain(SEXP);
extern SEXP xMain(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"clusteringMain", (DL_FUNC) &clusteringMain, 1},
    {"learnMain",      (DL_FUNC) &learnMain,      1},
    {"predictMain",    (DL_FUNC) &predictMain,    1},
    {"xMain",          (DL_FUNC) &xMain,          1},
    {NULL, NULL, 0}
};

void R_init_Rmixmod(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
