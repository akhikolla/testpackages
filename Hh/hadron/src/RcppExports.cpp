// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// read_nissa_textcf_kernel
NumericMatrix read_nissa_textcf_kernel(CharacterVector file_basenames_to_read, CharacterVector smear_combs_to_read, const unsigned int nts, DataFrame combs_to_read);
RcppExport SEXP _hadron_read_nissa_textcf_kernel(SEXP file_basenames_to_readSEXP, SEXP smear_combs_to_readSEXP, SEXP ntsSEXP, SEXP combs_to_readSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type file_basenames_to_read(file_basenames_to_readSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type smear_combs_to_read(smear_combs_to_readSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type nts(ntsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type combs_to_read(combs_to_readSEXP);
    rcpp_result_gen = Rcpp::wrap(read_nissa_textcf_kernel(file_basenames_to_read, smear_combs_to_read, nts, combs_to_read));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP alphas(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP cdh_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP cdhnew_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP invcosh(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_hadron_read_nissa_textcf_kernel", (DL_FUNC) &_hadron_read_nissa_textcf_kernel, 4},
    {"alphas",   (DL_FUNC) &alphas,    5},
    {"cdh_c",    (DL_FUNC) &cdh_c,    13},
    {"cdhnew_c", (DL_FUNC) &cdhnew_c, 11},
    {"invcosh",  (DL_FUNC) &invcosh,   5},
    {NULL, NULL, 0}
};

RcppExport void R_init_hadron(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
