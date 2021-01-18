// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fpa_optim
List fpa_optim(double N, double p, double beta, double eta, int maxiter, bool randEta, double gloMin, double objfit, double D, double Lower, double Upper, Function FUN);
RcppExport SEXP _MetaheuristicFPA_fpa_optim(SEXP NSEXP, SEXP pSEXP, SEXP betaSEXP, SEXP etaSEXP, SEXP maxiterSEXP, SEXP randEtaSEXP, SEXP gloMinSEXP, SEXP objfitSEXP, SEXP DSEXP, SEXP LowerSEXP, SEXP UpperSEXP, SEXP FUNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< bool >::type randEta(randEtaSEXP);
    Rcpp::traits::input_parameter< double >::type gloMin(gloMinSEXP);
    Rcpp::traits::input_parameter< double >::type objfit(objfitSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type Lower(LowerSEXP);
    Rcpp::traits::input_parameter< double >::type Upper(UpperSEXP);
    Rcpp::traits::input_parameter< Function >::type FUN(FUNSEXP);
    rcpp_result_gen = Rcpp::wrap(fpa_optim(N, p, beta, eta, maxiter, randEta, gloMin, objfit, D, Lower, Upper, FUN));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MetaheuristicFPA_fpa_optim", (DL_FUNC) &_MetaheuristicFPA_fpa_optim, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_MetaheuristicFPA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
