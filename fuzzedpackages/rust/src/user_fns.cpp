// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// User-supplied C++ functions for logf.

// Note that currently the only interface available in rust is
// double fun(const Rcpp::NumericVector& x, const Rcpp::List& pars).
// However, as shown in the function logdmvnorm below RcppArmadillo
// functions can be used inside the function.

// Each function must be prefaced by the line: // [[Rcpp::export]]

// One-dimensional standard normal.

// [[Rcpp::export]]
double logdN01(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return (-pow(x[0], 2.0) / 2.0) ;
}

// Two-dimensional normal with zero-mean and unit variances.

// [[Rcpp::export]]
double logdnorm2(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  double rho = pars["rho"] ;
  double div = 2.0 * (1.0 - pow(rho, 2.0)) ;
  return (-(pow(x[0], 2.0) - 2.0 * rho * x[0] * x[1] + pow(x[1], 2.0)) / div) ;
}

// d-dimensional normal with zero-mean and covariance matrix sigma.

// [[Rcpp::export]]
double logdmvnorm(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  arma::mat sigma = as<arma::mat>(pars["sigma"]) ;
  arma::vec y = Rcpp::as<arma::vec>(x) ;
  double qform = arma::as_scalar(y.t() * arma::inv(sigma) * y) ;
  return -qform / 2.0  ;
}

// Cauchy

// [[Rcpp::export]]
double logcauchy(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return -log(1 + pow(x[0], 2.0)) ;
}

// Half-Cauchy

// [[Rcpp::export]]
double loghalfcauchy(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  if (x[0] < 0)
    return R_NegInf ;
  return -log(1 + pow(x[0], 2.0)) ;
}

// Mixture of normals

// [[Rcpp::export]]
double lognormalmix(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  double mu = pars["mu"] ;
  double p = pars["p"] ;
  double comp1 = exp(-pow(x[0], 2.0) / 2.0)  ;
  double comp2 = exp(-pow(x[0] - mu, 2.0) / 2.0)  ;
  return log(p * comp1 + (1 - p) * comp2) ;
}

// Example 4 from Wakefield et al (1991).

// [[Rcpp::export]]
double lognormt(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  arma::vec mean = as<arma::vec>(pars["mean"]) ;
  arma::mat sigma1 = as<arma::mat>(pars["sigma1"]) ;
  arma::mat sigma2 = as<arma::mat>(pars["sigma2"]) ;
  arma::vec y1 = Rcpp::as<arma::vec>(x) ;
  arma::vec y = y1 - mean ;
  arma::vec y2 = Rcpp::as<arma::vec>(x) ;
  double log_h1 = -arma::as_scalar(y.t() * arma::inv(sigma1) * y) / 2 ;
  double log_h2 = -2 * log(1 + arma::as_scalar(y.t() * arma::inv(sigma1) * y)
                             / 2) ;
  return log_h1 + log_h2  ;
}

// Lognormal(mu, sigma).

// [[Rcpp::export]]
double logdlnorm(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  double mu = pars["mu"] ;
  double sigma = pars["sigma"] ;
  if (x[0] > 0)
    return -log(x[0]) - pow(log(x[0]) - mu, 2.0) / (2.0 * pow(sigma, 2.0)) ;
  else
    return R_NegInf ;
}

// Gamma(alpha, 1).

// [[Rcpp::export]]
double logdgamma(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  double shp = pars["alpha"] ;
  return Rcpp::dgamma(x, shp, 1.0, 1)[0] ;
}

// Generalized Pareto posterior based on an MDI prior truncated to
// shape parameter xi >= -1.

// [[Rcpp::export]]
double loggp(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  Rcpp::NumericVector gpd_data = ss["gpd_data"] ;
  int m = ss["m"] ;
  double xm = ss["xm"] ;
  double sum_gp = ss["sum_gp"] ;
  if (x[0] <= 0 || x[1] <= -x[0] / xm)
    return R_NegInf ;
  double loglik ;
  Rcpp::NumericVector sdat = gpd_data / x[0] ;
  Rcpp::NumericVector zz = 1 + x[1] * sdat ;
  if (std::abs(x[1]) > 1e-6) {
    loglik = -m * log(x[0]) - (1.0 + 1.0 / x[1]) * sum(log(zz)) ;
  } else {
    double t1, t2, sdatj ;
    double total = 0;
    for(int j = 0; j < m; ++j) {
      sdatj = sdat[j] ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        total += pow(-1.0, i) * t1 * t2 * pow(x[1], i) / i / (i + 1) ;
      }
    }
    loglik = -m * log(x[0]) - sum_gp / x[0] - total ;
  }
  // MDI prior.
  if (x[1] < -1)
    return R_NegInf ;
  double logprior = -log(x[0]) - x[1] - 1 ;
  return (logprior + loglik) ;
}

// A function to create external pointers to the functions to evaluate logf.
// See http://gallery.rcpp.org/articles/passing-cpp-function-pointers/
// If you write a new function above called new_name then add something
// like the following.
//
// else if (fstr == "new_name")
//   return(Rcpp::XPtr<funcPtr>(new funcPtr(&new_name))) ;

//' Create external pointer to a C++ function for logf
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP create_xptr(std::string fstr) {
  typedef double (*funcPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "logdN01")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&logdN01))) ;
  else if (fstr == "logdnorm2")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&logdnorm2))) ;
  else if (fstr == "logdmvnorm")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&logdmvnorm))) ;
  else if (fstr == "lognormt")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&lognormt))) ;
  else if (fstr == "loghalfcauchy")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&loghalfcauchy))) ;
  else if (fstr == "logcauchy")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&logcauchy))) ;
  else if (fstr == "lognormalmix")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&lognormalmix))) ;
  else if (fstr == "logdlnorm")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&logdlnorm))) ;
  else if (fstr == "logdgamma")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&logdgamma))) ;
  else if (fstr == "loggp")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&loggp))) ;
  else
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ;
}

// User-supplied C++ functions for log_j, the log-Jacobian of the
// transformation from theta to phi.

// [[Rcpp::export]]
double neglog(const Rcpp::NumericVector& theta,
              const Rcpp::List& user_args) {
  return -log(theta[0]) ;
}

// [[Rcpp::export]]
double bc_log_j(const Rcpp::NumericVector& theta,
                const Rcpp::List& user_args) {
  double lambda = user_args["lambda"] ;
  return (lambda - 1.0) * log(theta[0]) ;
}

// A function to create external pointers to functions to evaluate log_j.

//' Create external pointer to a C++ function for log_j
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP create_log_j_xptr(std::string fstr) {
  typedef double (*logjacPtr)(const Rcpp::NumericVector& theta,
                  const Rcpp::List& user_args) ;
  if (fstr == "neglog")
    return(Rcpp::XPtr<logjacPtr>(new logjacPtr(&neglog))) ;
  else if (fstr == "bc")
    return(Rcpp::XPtr<logjacPtr>(new logjacPtr(&bc_log_j))) ;
  else
    return(Rcpp::XPtr<logjacPtr>(R_NilValue)) ;
}

// User-supplied C++ functions for phi_to_theta, which performs
// transformation from phi to theta.

// [[Rcpp::export]]
Rcpp::NumericVector exptrans(const Rcpp::NumericVector& phi,
                                const Rcpp::List& user_args) {
  return exp(phi);
}

// [[Rcpp::export]]
Rcpp::NumericVector vecpower(const Rcpp::NumericVector& base,
                             const Rcpp::NumericVector& exp) {
  int n = base.size() ;
  Rcpp::NumericVector res(n) ;
  for(int i=0; i < n; i++) {
    res[i] = std::pow(base[i], exp[i]) ;
  }
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector bc_phi_to_theta(const Rcpp::NumericVector& phi,
                                    const Rcpp::List& user_args) {
  Rcpp::NumericVector lambda = user_args["lambda"] ;
  Rcpp::NumericVector test = phi * lambda + 1.0 ;
  return ifelse(test > 0, vecpower(test, 1.0 / lambda), NA_REAL) ;
}

// [[Rcpp::export]]
Rcpp::NumericVector gp_phi_to_theta(const Rcpp::NumericVector& phi,
                                    const Rcpp::List& user_args) {
  double xm = user_args["xm"] ;
  Rcpp::NumericVector val(2);
  val[0] = phi[0] ;
  val[1] = phi[1] - phi[0] / xm ;
  return val ;
}

// A function to create external pointers to functions to evaluate
// phi_to_theta.

//' Create external pointer to a C++ function for phi_to_theta
//'
//' @param fstr A string indicating the C++ function required.
//'
//' @export
// [[Rcpp::export]]
SEXP create_phi_to_theta_xptr(std::string fstr) {
  typedef Rcpp::NumericVector (*p2tPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  if (fstr == "exponential")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&exptrans))) ;
  else if (fstr == "bc")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&bc_phi_to_theta))) ;
  else if (fstr == "gp")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&gp_phi_to_theta))) ;
  else
    return(Rcpp::XPtr<p2tPtr>(R_NilValue)) ;
}

// We could create the external pointers when this file is sourced using
// this embedded R code below and/or (re)create them using the relevant
// pointer-creation functions in an R session or R package.

/*** R
ptr_N01 <- create_xptr("logdN01")
ptr_bvn <- create_xptr("logdnorm2")
ptr_mvn <- create_xptr("logdmvnorm")
ptr_c <- create_xptr("logcauchy")
ptr_hc <- create_xptr("loghalfcauchy")
ptr_normt <- create_xptr("lognormt")
ptr_lnorm <- create_xptr("logdlnorm")
ptr_gam <- create_xptr("logdgamma")
ptr_gp <- create_xptr("loggp")
ptr_phi_to_theta_lnorm <- create_phi_to_theta_xptr("exponential")
ptr_log_j_lnorm <- create_log_j_xptr("neglog")
ptr_phi_to_theta_bc <- create_phi_to_theta_xptr("bc")
ptr_log_j_bc <- create_log_j_xptr("bc")
ptr_phi_to_theta_gp <- create_phi_to_theta_xptr("gp")
*/

