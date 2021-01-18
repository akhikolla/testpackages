#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
Rcpp::List single_Fisher_test(arma::rowvec t,
                              bool correct, 
                              bool ret_all_probs ) {
  double n00 = t(0);
  double lo = t(10);
  vec support = linspace(lo, t(9), t(9)-lo+1);
  vec all_probs(support.n_elem);
  unsigned int c;
#pragma omp parallel for private(c)
  for (c=0; c<support.n_elem; c++) {
    all_probs(c) = R::dhyper(support(c), t(4), t(5), t(6), true);
  }

  all_probs = exp(all_probs - max(all_probs));
  all_probs /= sum(all_probs);
  double relErr = 1.0 + 1e-7;
  double pv = accu(all_probs(find(all_probs <= all_probs(n00-lo)*relErr)));
  if (pv==0) pv = DOUBLE_XMIN;
  double pv_correct = 0;
  if (correct) {
    double relErrNeg = 1 - 1e-7;
    pv_correct = (pv+accu(all_probs(find(all_probs <= all_probs(n00-lo)*relErrNeg))))/2;
    if (pv_correct==0) pv_correct = DOUBLE_XMIN;
  }
  if (correct & ret_all_probs) {
    all_probs = all_probs(find(all_probs>0));
    return Rcpp::List::create(
      Rcpp::Named( "all.probs" ) = all_probs,
      Rcpp::Named( "pv" ) = pv,
      Rcpp::Named( "pv.correct" ) = pv_correct  );
  }
  if (correct & !ret_all_probs) {
    return Rcpp::List::create(
      Rcpp::Named( "pv" ) = pv,
      Rcpp::Named( "pv.correct" ) = pv_correct  );
  }
  if (!correct & ret_all_probs) {
    all_probs = all_probs(find(all_probs>0));
    return Rcpp::List::create(
      Rcpp::Named( "all.probs" ) = all_probs,
      Rcpp::Named( "pv" ) = pv  );
  }
  if (!correct & !ret_all_probs) {
      return Rcpp::List::create( Rcpp::Named( "pv" ) = pv  );
  }
  
  return R_NilValue;
}
