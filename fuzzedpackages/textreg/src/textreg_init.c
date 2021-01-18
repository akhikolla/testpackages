#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP textreg_build_corpus(SEXP, SEXP, SEXP, SEXP);
extern SEXP textreg_textreg(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"textreg_build_corpus", (DL_FUNC) &textreg_build_corpus, 4},
    {"textreg_textreg",      (DL_FUNC) &textreg_textreg,      2},
    {NULL, NULL, 0}
};

void R_init_textreg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
