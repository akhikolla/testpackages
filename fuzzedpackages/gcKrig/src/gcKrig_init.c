#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP FHUBbinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP FHUBNB2(SEXP, SEXP, SEXP, SEXP);
extern SEXP FHUBNB2binomial(SEXP, SEXP, SEXP, SEXP);
extern SEXP FHUBZIP(SEXP, SEXP, SEXP, SEXP);
extern SEXP FHUBZIPbinomial(SEXP, SEXP, SEXP, SEXP);
extern SEXP FHUBZIPNB2(SEXP, SEXP, SEXP, SEXP);
extern SEXP ghkgcmr(SEXP, SEXP, SEXP, SEXP);
extern SEXP ghkgcmr2(SEXP, SEXP, SEXP, SEXP);
extern SEXP ldgc(SEXP, SEXP);
extern SEXP mvnintGHKcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mvnintGHKOcpp(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"FHUBbinom",       (DL_FUNC) &FHUBbinom,       4},
    {"FHUBNB2",         (DL_FUNC) &FHUBNB2,         4},
    {"FHUBNB2binomial", (DL_FUNC) &FHUBNB2binomial, 4},
    {"FHUBZIP",         (DL_FUNC) &FHUBZIP,         4},
    {"FHUBZIPbinomial", (DL_FUNC) &FHUBZIPbinomial, 4},
    {"FHUBZIPNB2",      (DL_FUNC) &FHUBZIPNB2,      4},
    {"ghkgcmr",         (DL_FUNC) &ghkgcmr,         4},
    {"ghkgcmr2",        (DL_FUNC) &ghkgcmr2,        4},
    {"ldgc",            (DL_FUNC) &ldgc,            2},
    {"mvnintGHKcpp",    (DL_FUNC) &mvnintGHKcpp,    5},
    {"mvnintGHKOcpp",   (DL_FUNC) &mvnintGHKOcpp,   5},
    {NULL, NULL, 0}
};

void R_init_gcKrig(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
