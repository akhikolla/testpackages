#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/


/* .Call calls */
extern SEXP CCK(SEXP, SEXP);
extern SEXP PVCk(SEXP, SEXP);
extern SEXP TIK(SEXP, SEXP, SEXP, SEXP);
extern SEXP TIKlc(SEXP, SEXP, SEXP, SEXP);
extern SEXP TIK1(SEXP, SEXP, SEXP);
extern SEXP TSAVE(SEXP, SEXP, SEXP, SEXP);
extern SEXP TSIR(SEXP, SEXP, SEXP, SEXP);
extern SEXP lblinM(SEXP, SEXP);
extern SEXP lbsqM(SEXP, SEXP);
extern SEXP EIGEN(SEXP);
extern SEXP PREPBSS(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"CCK",    (DL_FUNC) &CCK,    2},
    {"PVCk",   (DL_FUNC) &PVCk,   2},
    {"TIK",    (DL_FUNC) &TIK,    4},
    {"TIKlc",  (DL_FUNC) &TIKlc,  4},
    {"TIK1",   (DL_FUNC) &TIK1,   3},
    {"TSAVE",  (DL_FUNC) &TSAVE,  4},
    {"TSIR",   (DL_FUNC) &TSIR,   4},
    {"lblinM", (DL_FUNC) &lblinM, 2},
    {"lbsqM",  (DL_FUNC) &lbsqM,  2},
    {"EIGEN",  (DL_FUNC) &EIGEN,  1},
    {"PREPBSS",  (DL_FUNC) &PREPBSS,  2},
    {NULL, NULL, 0}
};



void R_init_tsBSS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
