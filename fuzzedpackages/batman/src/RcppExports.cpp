// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// to_logical_
LogicalVector to_logical_(std::vector < std::string > x, std::vector < std::string > trues, std::vector < std::string > falses);
RcppExport SEXP _batman_to_logical_(SEXP xSEXP, SEXP truesSEXP, SEXP falsesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector < std::string > >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::vector < std::string > >::type trues(truesSEXP);
    Rcpp::traits::input_parameter< std::vector < std::string > >::type falses(falsesSEXP);
    rcpp_result_gen = Rcpp::wrap(to_logical_(x, trues, falses));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_batman_to_logical_", (DL_FUNC) &_batman_to_logical_, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_batman(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
