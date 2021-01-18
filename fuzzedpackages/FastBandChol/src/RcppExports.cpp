#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// BandCholcpp
List BandCholcpp(arma::mat x, int k);
RcppExport SEXP FastBandChol_BandCholcpp(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    __result = Rcpp::wrap(BandCholcpp(x, k));
    return __result;
END_RCPP
}
