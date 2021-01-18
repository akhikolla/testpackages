#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP DDRTree_DDRTree_reduce_dim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP DDRTree_pca_projection(SEXP, SEXP);
extern SEXP DDRTree_sqdist(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"DDRTree_DDRTree_reduce_dim", (DL_FUNC) &DDRTree_DDRTree_reduce_dim, 12},
    {"DDRTree_pca_projection",     (DL_FUNC) &DDRTree_pca_projection,      2},
    {"DDRTree_sqdist",             (DL_FUNC) &DDRTree_sqdist,              2},
    {NULL, NULL, 0}
};

void R_init_DDRTree(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
