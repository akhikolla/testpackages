// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kde_estimate
NumericVector kde_estimate(NumericMatrix fishnet, NumericMatrix points, double bw, String kernel, bool scaled, double decay, NumericVector weights);
RcppExport SEXP _SpatialKDE_kde_estimate(SEXP fishnetSEXP, SEXP pointsSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP scaledSEXP, SEXP decaySEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type fishnet(fishnetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< String >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< bool >::type scaled(scaledSEXP);
    Rcpp::traits::input_parameter< double >::type decay(decaySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(kde_estimate(fishnet, points, bw, kernel, scaled, decay, weights));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SpatialKDE_kde_estimate", (DL_FUNC) &_SpatialKDE_kde_estimate, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SpatialKDE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
