#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP hyper2_accessor(SEXP, SEXP, SEXP);
extern SEXP hyper2_addL(SEXP, SEXP, SEXP, SEXP);
extern SEXP hyper2_assigner(SEXP, SEXP, SEXP, SEXP);
extern SEXP hyper2_differentiate(SEXP, SEXP, SEXP, SEXP);
extern SEXP hyper2_equality(SEXP, SEXP, SEXP, SEXP);
extern SEXP hyper2_evaluate(SEXP, SEXP, SEXP);
extern SEXP hyper2_identityL(SEXP, SEXP);
extern SEXP hyper2_overwrite(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"hyper2_accessor",      (DL_FUNC) &hyper2_accessor,      3},
    {"hyper2_addL",          (DL_FUNC) &hyper2_addL,          4},
    {"hyper2_assigner",      (DL_FUNC) &hyper2_assigner,      4},
    {"hyper2_differentiate", (DL_FUNC) &hyper2_differentiate, 4},
    {"hyper2_equality",      (DL_FUNC) &hyper2_equality,      4},
    {"hyper2_evaluate",      (DL_FUNC) &hyper2_evaluate,      3},
    {"hyper2_identityL",     (DL_FUNC) &hyper2_identityL,     2},
    {"hyper2_overwrite",     (DL_FUNC) &hyper2_overwrite,     4},
    {NULL, NULL, 0}
};

void R_init_hyper2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
