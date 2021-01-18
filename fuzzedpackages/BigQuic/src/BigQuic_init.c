#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP BigQuic_BigQuicHelper(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"BigQuic_BigQuicHelper", (DL_FUNC) &BigQuic_BigQuicHelper, 1},
    {NULL, NULL, 0}
};

void R_init_BigQuic(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
