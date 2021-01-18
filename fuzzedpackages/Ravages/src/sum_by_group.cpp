#include <Rcpp.h>
#include <iostream>
#include <ctime>

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector sum_by_group(NumericVector stat, IntegerVector group) {
  int nlevels = as<CharacterVector>(group.attr("levels")).size();
  int n = stat.size();
  if(group.size() != n) {
    stop("stat and group don't have the same length");
  }
  NumericVector R(nlevels);
  for(int i = 0; i < n; i++) 
    if(!NumericVector::is_na(stat[i]) && !IntegerVector::is_na(group[i])) R[ group[i]-1 ] += stat[i];
  
  return R;
}


RcppExport SEXP oz_sum_by_group(SEXP statSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_by_group(stat, group));
    return rcpp_result_gen;
END_RCPP
}

