#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void locLastZero(void *, void *, void *);
extern void nocopy_kmc_data(void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP kmcRCPP_KMCDATA(SEXP, SEXP, SEXP, SEXP);
extern SEXP kmcRCPP_RevCHECK(SEXP);
extern SEXP kmcomegalambda(SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"locLastZero",     (DL_FUNC) &locLastZero,     3},
    {"nocopy_kmc_data", (DL_FUNC) &nocopy_kmc_data, 5},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"kmcRCPP_KMCDATA",  (DL_FUNC) &kmcRCPP_KMCDATA,  4},
    {"kmcRCPP_RevCHECK", (DL_FUNC) &kmcRCPP_RevCHECK, 1},
    {"kmcomegalambda",   (DL_FUNC) &kmcomegalambda,   4},
    {NULL, NULL, 0}
};

void R_init_kmc(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
