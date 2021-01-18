// generated using tools::package_native_routine_registration_skeleton(".", character_only = FALSE)

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _sboost_adaboost(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _sboost_get_contingency_cpp(SEXP, SEXP, SEXP);
extern SEXP _sboost_predict_cpp(SEXP, SEXP);
extern SEXP _sboost_score_classifier_features_cpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_sboost_adaboost",                      (DL_FUNC) &_sboost_adaboost,                      6},
    {"_sboost_get_contingency_cpp",           (DL_FUNC) &_sboost_get_contingency_cpp,           3},
    {"_sboost_predict_cpp",                   (DL_FUNC) &_sboost_predict_cpp,                   2},
    {"_sboost_score_classifier_features_cpp", (DL_FUNC) &_sboost_score_classifier_features_cpp, 2},
    {NULL, NULL, 0}
};

void R_init_sboost(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
