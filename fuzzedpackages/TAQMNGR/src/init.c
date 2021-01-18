	
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP Aggregate(SEXP, SEXP, SEXP, SEXP);
extern SEXP CleaningReport(SEXP, SEXP);
extern SEXP CleanTickByTick(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Import(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"Aggregate",       (DL_FUNC) &Aggregate,        4},
    {"CleaningReport",  (DL_FUNC) &CleaningReport,   2},
    {"CleanTickByTick", (DL_FUNC) &CleanTickByTick,  6},
    {"Import",          (DL_FUNC) &Import,          10},
    {NULL, NULL, 0}
};

void R_init_TAQMNGR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
