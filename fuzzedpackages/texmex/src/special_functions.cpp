#include <Rcpp.h>

#include "special_functions.h"


//' Accurately compute (exp(x) - 1) / x
//' @param x numeric vector
//' @return numeric vector
// [[Rcpp::export(name='.exprel', rng=FALSE)]]
Rcpp::NumericVector warp_dexprl(const Rcpp::NumericVector& x) {
  return Rcpp::sapply(x, dexprl);
}

//' Accurately compute log(1 + x) / x
//' @param x numeric vector
//' @return numeric vector
// [[Rcpp::export(name='.log1prel', rng=FALSE)]]
Rcpp::NumericVector wrap_log1prel(const Rcpp::NumericVector& x) {
  return Rcpp::sapply(x, log1prel);
}

//' Accurately compute log(1-exp(x))
//' @param x numeric vector
//' @return a numeric vector
// [[Rcpp::export(name=".log1mexp", rng=FALSE)]]
Rcpp::NumericVector wrap_log1mexp(const Rcpp::NumericVector& x) {
  return Rcpp::sapply(x, texmex_log1mexp);
}

//' Compute pmax(x y, -1) in such a way that zeros in x beat
//' infinities in y.
//'
//' This is a common pattern in much of the distribution code, so it's
//' worth factoring out.
//' @param x a numeric vector
//' @param y a numeric vector
//' @return an appropriate numeric vector
// [[Rcpp::export(name=".specfun.safe.product", rng=FALSE)]]
Rcpp::NumericVector wrap_safe_product(const Rcpp::NumericVector &x,
				      const Rcpp::NumericVector &y) {
  const R_xlen_t size = std::max(x.size(), y.size());
  return Rcpp::mapply(Rcpp::rep_len(x, size),
		      Rcpp::rep_len(y, size),
		      safe_product);
}
