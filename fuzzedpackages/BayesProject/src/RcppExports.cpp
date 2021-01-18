// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// bayes_vhat
MatrixXd bayes_vhat(MatrixXd x, VectorXd timePoints, double K);
RcppExport SEXP _BayesProject_bayes_vhat(SEXP xSEXP, SEXP timePointsSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type timePoints(timePointsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(bayes_vhat(x, timePoints, K));
    return rcpp_result_gen;
END_RCPP
}
// bayes_cpt
VectorXd bayes_cpt(MatrixXd x, VectorXd timePoints, double K);
RcppExport SEXP _BayesProject_bayes_cpt(SEXP xSEXP, SEXP timePointsSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type timePoints(timePointsSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(bayes_cpt(x, timePoints, K));
    return rcpp_result_gen;
END_RCPP
}
// sum_max_cusum
VectorXd sum_max_cusum(MatrixXd x, bool sum_cusum);
RcppExport SEXP _BayesProject_sum_max_cusum(SEXP xSEXP, SEXP sum_cusumSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type sum_cusum(sum_cusumSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_max_cusum(x, sum_cusum));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesProject_bayes_vhat", (DL_FUNC) &_BayesProject_bayes_vhat, 3},
    {"_BayesProject_bayes_cpt", (DL_FUNC) &_BayesProject_bayes_cpt, 3},
    {"_BayesProject_sum_max_cusum", (DL_FUNC) &_BayesProject_sum_max_cusum, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesProject(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
