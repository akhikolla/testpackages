
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP fastJT_fastJT(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fastJT_fastJTmp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"fastJT_fastJT",   (DL_FUNC) &fastJT_fastJT,   5},
    {"fastJT_fastJTmp", (DL_FUNC) &fastJT_fastJTmp, 6},
    {NULL, NULL, 0}
};

void R_init_fastJT(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

