#include "partial_sort_visitor.h"

// [[Rcpp::export]]
SEXP partial_sort_cpp(SEXP rv, SEXP rpivots) {
  Rcpp::XPtr<ldat::vec> v(rv);
  // convert R vector of pivots to std::vector
  Rcpp::NumericVector pivots_r(rpivots);
  std::vector<ldat::vec::vecsize> pivots;
  for (R_xlen_t i = 0; i < pivots_r.length(); ++i) {
    if (ldat::is_na(pivots_r[i]))
      throw Rcpp::exception("Missing values in the pivots");
    pivots.push_back(pivots_r[i] - 1);
  }
  // call visitor
  psort_visitor visitor(pivots);
  v->visit(&visitor);
  return R_NilValue;
}

