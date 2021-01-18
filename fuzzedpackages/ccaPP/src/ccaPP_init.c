#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP R_ccaPP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_corKendall(SEXP, SEXP, SEXP);
extern SEXP R_corM(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_corPearson(SEXP, SEXP);
extern SEXP R_corQuadrant(SEXP, SEXP, SEXP);
extern SEXP R_corSpearman(SEXP, SEXP, SEXP);
extern SEXP R_fastMAD(SEXP, SEXP);
extern SEXP R_fastMedian(SEXP);
extern SEXP R_l1Median(SEXP);
extern SEXP R_rank(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_ccaPP",       (DL_FUNC) &R_ccaPP,       9},
    {"R_corKendall",  (DL_FUNC) &R_corKendall,  3},
    {"R_corM",        (DL_FUNC) &R_corM,        5},
    {"R_corPearson",  (DL_FUNC) &R_corPearson,  2},
    {"R_corQuadrant", (DL_FUNC) &R_corQuadrant, 3},
    {"R_corSpearman", (DL_FUNC) &R_corSpearman, 3},
    {"R_fastMAD",     (DL_FUNC) &R_fastMAD,     2},
    {"R_fastMedian",  (DL_FUNC) &R_fastMedian,  1},
    {"R_l1Median",    (DL_FUNC) &R_l1Median,    1},
    {"R_rank",        (DL_FUNC) &R_rank,        1},
    {NULL, NULL, 0}
};

void R_init_ccaPP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
