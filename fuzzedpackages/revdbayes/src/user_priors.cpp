// [[Rcpp::depends(Rcpp)]]

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// Generalized Pareto log-priors

// [[Rcpp::export]]
double user_gp_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi)
    return R_NegInf ;
  return -log(x[0]) ;
}

// Generalized Extreme Value log-priors

// [[Rcpp::export]]
double user_gev_norm(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  if (x[1] <= 0)
    return R_NegInf ;
  Rcpp::NumericVector mean = ppars["mean"] ;
  Rcpp::NumericVector icov = ppars["icov"] ;
  double c0 = x[0] - mean[0] ;
  double c1 = log(x[1]) - mean[1] ;
  double c2 = x[2] - mean[2] ;
  double ld = icov[0] * pow(c0, 2.0) + 2 * icov[1] * c0 * c1 +
    2 * icov[2] * c0 * c2 + icov[3] * pow(c1, 2.0) + 2 * icov[4] * c1 * c2 +
    icov[5] * pow(c2, 2.0) ;
  return (-ld / 2 - log(x[1])) ;
}

// [[Rcpp::export]]
double user_gev_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  return -log(x[1]) ;
}

//' Create an external pointer to a C++ prior
//'
//' This function provides an example of a way in which a user can specify
//' their own prior density to \code{\link{rpost_rcpp}}.
//' More specifically, a function like this (the user will need to create
//' an edited version tailored to their own C++ function(s)) can be used to
//' generate an external pointer to a compiled C++ function that evaluates
//' the log-prior density.  Please see the vignette
//' "Faster simulation using revdbayes" for more information.
//'
//' @details
//' Suppose that the user's C++ functions are in a file called "user_fns.cpp".
//' These functions must be compiled and made available to R before the
//' pointer is created. This can be achieved using the function
//' \code{\link[Rcpp]{sourceCpp}} in the \strong{Rcpp} package
//' or using RStudio's Source button on the editor toolbar.
//'
//' For details see the examples in the documentation of the functions
//' \code{\link{rpost_rcpp}} and \code{\link{set_prior}},
//' the vignette "Faster simulation using revdbayes"
//' and the vignette "Rusting Faster: Simulation using Rcpp" in the package
//' \strong{rust}.
//'
//' @param fstr A string indicating the C++ function required.
//' @return An external pointer.
//' @seealso \code{\link{set_prior}} to specify a prior distribution using
//'   an external pointer returned by \code{create_prior_xptr} and for
//'   details of in-built named prior distributions.
//' @seealso The examples in the documentation of \code{\link{rpost_rcpp}}.
//' @examples
//' ptr_gp_flat <- create_prior_xptr("gp_flat")
//' prior_cfn <- set_prior(prior = ptr_gp_flat, model = "gp", min_xi = -1)
//'
//' ptr_gev_flat <- create_prior_xptr("gev_flat")
//' prior_cfn <- set_prior(prior = ptr_gev_flat, model = "gev", min_xi = -1,
//'                        max_xi = Inf)
//'
//' mat <- diag(c(10000, 10000, 100))
//' ptr_gev_norm <- create_prior_xptr("gev_norm")
//' pn_u <- set_prior(prior = ptr_gev_norm, model = "gev", mean = c(0,0,0),
//'                   icov = solve(mat))
//' @export
// [[Rcpp::export]]
SEXP create_prior_xptr(std::string fstr) {
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  if (fstr == "gp_flat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&user_gp_flat))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&user_gev_norm))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<priorPtr>(new priorPtr(&user_gev_flat))) ;
  else
    return(Rcpp::XPtr<priorPtr>(R_NilValue)) ;
}

// We could create the external pointers when this file is sourced using
// this embedded R code below and/or (re)create them using the relevant
// pointer-creation functions in an R session or R package.

/*** R
  ptr_gp_flat <- create_prior_xptr("gp_flat")
  ptr_gev_norm <- create_prior_xptr("gev_norm")
  ptr_gev_flat <- create_prior_xptr("gev_flat")
*/
