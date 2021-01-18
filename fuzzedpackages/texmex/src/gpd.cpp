#include <algorithm>

#include <Rcpp.h>

#include "exponentiate.h"
#include "mapply.h"
#include "special_functions.h"


namespace {
  struct dgpd_func {
    typedef double result_type;

    inline double
    operator()(const double x,
	       const double sigma,
	       const double xi,
	       const double u) const {
      if (x < u) {
	return R_NegInf;
      }
      const double x_use = std::max((x - u) / sigma, 0.0);
      const double xi_x  = safe_product(xi, x_use);
      if (xi_x == -1.0) {
	return R_NegInf;
      }
      
      return -std::log(sigma) - log1p(xi_x) - log1prel(xi_x) * x_use;
    };
  };
}

// [[Rcpp::export(name=".dgpd", rng=FALSE)]]
Rcpp::NumericVector wrap_dgpd(const Rcpp::NumericVector &x,
			      const Rcpp::NumericVector &sigma,
			      const Rcpp::NumericVector &xi,
			      const Rcpp::NumericVector &u,
			      const bool log_d=false) {

  const R_xlen_t size = std::max(std::max(x.size(), sigma.size()),
				 std::max(xi.size(), u.size()));
  
  return perhaps_exp(mapply(Rcpp::rep_len(x, size),
			    Rcpp::rep_len(sigma, size),
			    Rcpp::rep_len(xi, size),
			    Rcpp::rep_len(u, size),
			    dgpd_func()),
		     log_d);
}
