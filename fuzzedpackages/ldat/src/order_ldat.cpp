#include "compare_visitor.h"

// [[Rcpp::export]]
SEXP order_ldat_cpp(Rcpp::List veclist) {
  // create vector with lvecs on which to sort
  std::vector<ldat::vec*> vecs;
  for (int i = 0; i < veclist.length(); ++i) {
    vecs.push_back(
      Rcpp::XPtr<ldat::vec>(static_cast<SEXP>(veclist[i]))
    );
  }
  if (vecs.size() < 1) 
    throw Rcpp::exception("No columns which to order.");
  // initialise result with current order
  ldat::vec::vecsize size = vecs[0]->size();
  std::unique_ptr<ldat::lvec<double> > result(new ldat::lvec<double>(size));
  for (ldat::vec::vecsize i = 0; i < size; ++i) result->set(i, i+1);
  // we are now going to sort the current order; but when comparing
  // we use the values pointed to by the order vector
  std::sort(result->begin(), result->end(), ldat_compare(vecs));
  return Rcpp::XPtr<ldat::vec>(result.release(), true);
}

