/* NOTE: this code is generated with tools::package_native_routine_registration_skeleton("R") */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _PRIMME_dprimme_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PRIMME_dprimme_svds_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PRIMME_primme_free_rcpp(SEXP);
extern SEXP _PRIMME_primme_get_member_rcpp(SEXP, SEXP);
extern SEXP _PRIMME_primme_initialize_rcpp();
extern SEXP _PRIMME_primme_set_member_rcpp(SEXP, SEXP, SEXP);
extern SEXP _PRIMME_primme_set_method_rcpp(SEXP, SEXP);
extern SEXP _PRIMME_primme_svds_free_rcpp(SEXP);
extern SEXP _PRIMME_primme_svds_get_member_rcpp(SEXP, SEXP);
extern SEXP _PRIMME_primme_svds_initialize_rcpp();
extern SEXP _PRIMME_primme_svds_set_member_rcpp(SEXP, SEXP, SEXP);
extern SEXP _PRIMME_primme_svds_set_method_rcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _PRIMME_zprimme_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _PRIMME_zprimme_svds_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_PRIMME_dprimme_rcpp",                (DL_FUNC) &_PRIMME_dprimme_rcpp,                7},
    {"_PRIMME_dprimme_svds_rcpp",           (DL_FUNC) &_PRIMME_dprimme_svds_rcpp,           7},
    {"_PRIMME_primme_free_rcpp",            (DL_FUNC) &_PRIMME_primme_free_rcpp,            1},
    {"_PRIMME_primme_get_member_rcpp",      (DL_FUNC) &_PRIMME_primme_get_member_rcpp,      2},
    {"_PRIMME_primme_initialize_rcpp",      (DL_FUNC) &_PRIMME_primme_initialize_rcpp,      0},
    {"_PRIMME_primme_set_member_rcpp",      (DL_FUNC) &_PRIMME_primme_set_member_rcpp,      3},
    {"_PRIMME_primme_set_method_rcpp",      (DL_FUNC) &_PRIMME_primme_set_method_rcpp,      2},
    {"_PRIMME_primme_svds_free_rcpp",       (DL_FUNC) &_PRIMME_primme_svds_free_rcpp,       1},
    {"_PRIMME_primme_svds_get_member_rcpp", (DL_FUNC) &_PRIMME_primme_svds_get_member_rcpp, 2},
    {"_PRIMME_primme_svds_initialize_rcpp", (DL_FUNC) &_PRIMME_primme_svds_initialize_rcpp, 0},
    {"_PRIMME_primme_svds_set_member_rcpp", (DL_FUNC) &_PRIMME_primme_svds_set_member_rcpp, 3},
    {"_PRIMME_primme_svds_set_method_rcpp", (DL_FUNC) &_PRIMME_primme_svds_set_method_rcpp, 4},
    {"_PRIMME_zprimme_rcpp",                (DL_FUNC) &_PRIMME_zprimme_rcpp,                7},
    {"_PRIMME_zprimme_svds_rcpp",           (DL_FUNC) &_PRIMME_zprimme_svds_rcpp,           7},
    {NULL, NULL, 0}
};

void R_init_PRIMME(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
