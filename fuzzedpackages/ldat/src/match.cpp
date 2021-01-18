#include "match_visitor.h"

// [[Rcpp::export]]
SEXP lmatch_cpp(SEXP rv, SEXP rvo, SEXP rtab, SEXP rtabo, SEXP rna_incomp) {
  Rcpp::XPtr<ldat::vec> v(rv);
  Rcpp::XPtr<ldat::vec> vo(rvo);
  if (v->size() != vo->size()) 
    throw Rcpp::exception("Lengths of vector and order of vector are unequal.");
  // check and convert table
  Rcpp::XPtr<ldat::vec> tab(rtab);
  Rcpp::XPtr<ldat::vec> tabo(rtabo);
  if (tab->size() != tabo->size()) 
    throw Rcpp::exception("Lengths of table and order of table are unequal.");
  bool na_incomp = Rcpp::as<bool>(rna_incomp);
  // call visitor
  match_visitor visitor(vo, tab, tabo, na_incomp);
  v->visit(&visitor);
  return Rcpp::XPtr<ldat::vec>(visitor.result(), true);
}

