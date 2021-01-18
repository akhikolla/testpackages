#include <Rcpp.h>

using namespace Rcpp;

List	omegalambda(SEXP kmctime,SEXP delta,SEXP lambda,SEXP gtmat);

RcppExport SEXP kmcomegalambda(SEXP kmctime,SEXP delta,SEXP lambda,SEXP gtmat) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = omegalambda( kmctime,delta,lambda,gtmat);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}


List RCPP_KMCDATA(SEXP kmctime,SEXP delta,SEXP lambda,SEXP gtmat);

RcppExport SEXP kmcRCPP_KMCDATA(SEXP kmctime,SEXP delta,SEXP lambda,SEXP gtmat) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = RCPP_KMCDATA( kmctime,delta,lambda,gtmat);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

List RevCHECK(SEXP x);

RcppExport SEXP kmcRCPP_RevCHECK(SEXP gmat){
BEGIN_RCPP
    SEXP __sexp_result;
  {
    List __result = RevCHECK(gmat);
    PROTECT(__sexp_result = Rcpp::wrap(__result));
  }
  UNPROTECT(1);
  return __sexp_result;
END_RCPP
}
