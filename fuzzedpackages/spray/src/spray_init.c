#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP spray_spray_accessor(SEXP, SEXP, SEXP);
extern SEXP spray_spray_add(SEXP, SEXP, SEXP, SEXP);
extern SEXP spray_spray_asum_exclude(SEXP, SEXP, SEXP);
extern SEXP spray_spray_asum_include(SEXP, SEXP, SEXP);
extern SEXP spray_spray_deriv(SEXP, SEXP, SEXP);
extern SEXP spray_spray_equality(SEXP, SEXP, SEXP, SEXP);
extern SEXP spray_spray_maker(SEXP, SEXP);
extern SEXP spray_spray_mult(SEXP, SEXP, SEXP, SEXP);
extern SEXP spray_spray_overwrite(SEXP, SEXP, SEXP, SEXP);
extern SEXP spray_spray_pmax(SEXP, SEXP, SEXP, SEXP);
extern SEXP spray_spray_pmin(SEXP, SEXP, SEXP, SEXP);
extern SEXP spray_spray_power(SEXP, SEXP, SEXP);
extern SEXP spray_spray_setter(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"spray_spray_accessor",     (DL_FUNC) &spray_spray_accessor,     3},
    {"spray_spray_add",          (DL_FUNC) &spray_spray_add,          4},
    {"spray_spray_asum_exclude", (DL_FUNC) &spray_spray_asum_exclude, 3},
    {"spray_spray_asum_include", (DL_FUNC) &spray_spray_asum_include, 3},
    {"spray_spray_deriv",        (DL_FUNC) &spray_spray_deriv,        3},
    {"spray_spray_equality",     (DL_FUNC) &spray_spray_equality,     4},
    {"spray_spray_maker",        (DL_FUNC) &spray_spray_maker,        2},
    {"spray_spray_mult",         (DL_FUNC) &spray_spray_mult,         4},
    {"spray_spray_overwrite",    (DL_FUNC) &spray_spray_overwrite,    4},
    {"spray_spray_pmax",         (DL_FUNC) &spray_spray_pmax,         4},
    {"spray_spray_pmin",         (DL_FUNC) &spray_spray_pmin,         4},
    {"spray_spray_power",        (DL_FUNC) &spray_spray_power,        3},
    {"spray_spray_setter",       (DL_FUNC) &spray_spray_setter,       4},
    {NULL, NULL, 0}
};

void R_init_spray(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
