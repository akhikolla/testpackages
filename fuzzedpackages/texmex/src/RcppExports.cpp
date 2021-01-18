// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// wrap_dgpd
Rcpp::NumericVector wrap_dgpd(const Rcpp::NumericVector& x, const Rcpp::NumericVector& sigma, const Rcpp::NumericVector& xi, const Rcpp::NumericVector& u, const bool log_d);
RcppExport SEXP _texmex_wrap_dgpd(SEXP xSEXP, SEXP sigmaSEXP, SEXP xiSEXP, SEXP uSEXP, SEXP log_dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type u(uSEXP);
    Rcpp::traits::input_parameter< const bool >::type log_d(log_dSEXP);
    rcpp_result_gen = Rcpp::wrap(wrap_dgpd(x, sigma, xi, u, log_d));
    return rcpp_result_gen;
END_RCPP
}
// warp_dexprl
Rcpp::NumericVector warp_dexprl(const Rcpp::NumericVector& x);
RcppExport SEXP _texmex_warp_dexprl(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(warp_dexprl(x));
    return rcpp_result_gen;
END_RCPP
}
// wrap_log1prel
Rcpp::NumericVector wrap_log1prel(const Rcpp::NumericVector& x);
RcppExport SEXP _texmex_wrap_log1prel(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(wrap_log1prel(x));
    return rcpp_result_gen;
END_RCPP
}
// wrap_log1mexp
Rcpp::NumericVector wrap_log1mexp(const Rcpp::NumericVector& x);
RcppExport SEXP _texmex_wrap_log1mexp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(wrap_log1mexp(x));
    return rcpp_result_gen;
END_RCPP
}
// wrap_safe_product
Rcpp::NumericVector wrap_safe_product(const Rcpp::NumericVector& x, const Rcpp::NumericVector& y);
RcppExport SEXP _texmex_wrap_safe_product(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(wrap_safe_product(x, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_texmex_wrap_dgpd", (DL_FUNC) &_texmex_wrap_dgpd, 5},
    {"_texmex_warp_dexprl", (DL_FUNC) &_texmex_warp_dexprl, 1},
    {"_texmex_wrap_log1prel", (DL_FUNC) &_texmex_wrap_log1prel, 1},
    {"_texmex_wrap_log1mexp", (DL_FUNC) &_texmex_wrap_log1mexp, 1},
    {"_texmex_wrap_safe_product", (DL_FUNC) &_texmex_wrap_safe_product, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_texmex(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
