// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ced_enc_detect
SEXP ced_enc_detect(SEXP x, SEXP enc_hint, SEXP lang_hint);
RcppExport SEXP _ced_ced_enc_detect(SEXP xSEXP, SEXP enc_hintSEXP, SEXP lang_hintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type enc_hint(enc_hintSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lang_hint(lang_hintSEXP);
    rcpp_result_gen = Rcpp::wrap(ced_enc_detect(x, enc_hint, lang_hint));
    return rcpp_result_gen;
END_RCPP
}
// ced_version
Rcpp::List ced_version();
RcppExport SEXP _ced_ced_version() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    rcpp_result_gen = Rcpp::wrap(ced_version());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ced_ced_enc_detect", (DL_FUNC) &_ced_ced_enc_detect, 3},
    {"_ced_ced_version", (DL_FUNC) &_ced_ced_version, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_ced(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
