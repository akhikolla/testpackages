#include <iostream>
#include "RcppArmadillo.h"

using namespace arma;
using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// =============================================================================
// ------------------------------ survivlas ------------------------------------
// =============================================================================
// Weibull survival function

// Weibull survival function in proportional hazard parametrisation
// The parametrisation assumed here is the Proportional hazard
// form used in McShane(2008)
//  \deqn{F(t) = exp(-\lambda * t^{\beta}}.
// where \eqn{\lambda} is the scale parameter and \eqn{\beta} is the shape
// parameter
// @param t double time point where the survival is to be evaluated at.
// @param distPars Rcpp::List with slots \code{scale} (the scale) and
// \code{shape} (the shape) parameters.
// @return a double giving the value of the survival function at time point
// \code{t} at the parameters value.
//
//' @keywords internal
// [[Rcpp::export]]
double sWeibull(double t, const Rcpp::List distPars) {
  double lambda = as <double> (distPars["scale"]);
  double beta = as <double> (distPars["shape"]);

  return exp( -lambda * pow( t, beta));
}

// Burr survival function

// Burr type X|| survival function in proportional hazard parametrisation
// The parametrisation assumed here is the Proportional hazard
// form
//  \deqn{F(t) = \frac{1}{(1 + \lambda t^\beta)^{\nu}}
// where \eqn{\lambda} is the scale parameter and \eqn{\beta} and \eqn{\nu}
// are the shape parameters.
// @param t double time point where the survival is to be evaluated at.
// @param distPars Rcpp::List with slots \code{scale} (the scale) and
// \code{shape1} (the shape1) and \code{shape2} (the shape2)parameters.
// @return a double giving the value of the survival function at time point
// \code{t} at the parameters value.
//
//' @keywords internal
// [[Rcpp::export]]
double sBurr(double t, const Rcpp::List distPars) {
  double lambda = as <double> (distPars["scale"]);
  double beta = as <double> (distPars["shape1"]);
  double nu = as <double> (distPars["shape2"]);

  double temp = 1 + lambda * pow(t, beta);
  return(1.0 / pow(temp, nu));
}

// gamma survival function

// gamma survival function in shape scale parametrisation.
//
// The \code{R::pgamma} is used to generate the survival function
// @param t double time point where the survival is to be evaluated at.
// @param distPars Rcpp::List with slots \code{scale} (the scale) and
// \code{shape} (the shape) parameter.
// @return a double giving the value of the survival function at time point
// \code{t} at the parameters value.
//
//' @keywords internal
// [[Rcpp::export]]
double sgamma(double t, const Rcpp::List distPars) {
  double shape = as <double> (distPars["shape"]);
  double rate  = as <double> (distPars["rate"]);

  return R::pgamma(t, shape, 1.0 / rate, false, false);
}

// generalized gamma survival function

// generalized gamma survival function in the Prentice parametrization.
//
// The Prentice parametrization is used which is preferred to the original
// parameterisation of the generalized gamma by Stacy (1962) since it is more
// numerically stable near to \eqn{Q = 0} (the log-normal distribution), and
// allows \eqn{Q \leq 0}.
// @param t double time point where the survival is to be evaluated at.
// @param distPars Rcpp::List with slots \code{mu} (location), \code{sigma}
// (scale parameters} and \code{Q} (the shape) parameter.
// @return a double giving the value of the survival function at time point
// \code{t} at the parameters value.
//
//' @keywords internal
// [[Rcpp::export]]
double sgengamma(double t, const Rcpp::List distPars) {
  double mu = as <double> (distPars["mu"]);
  double sigma = as <double> (distPars["sigma"]);
  double Q = as <double> (distPars["Q"]);

  double res = (log(t) + mu) / sigma; // modified to match glm-Poisson
  double Q2 = pow(Q, 2);
  double invQ2 = 1.0 / Q2;
  double tempGam, out;

  if (Q == 0.0) {
    out = 1 - R::pnorm(res, 0, 1, true, false);
  } else {
    tempGam = R::pgamma(exp(Q * res) * invQ2, invQ2, 1, true, false);
    if (Q > 0) {
      out = 1 - tempGam;
    } else {
      out = tempGam;
    }
  }
  return out ;
}

//' Wrapper to built in survival functions
//'
//' Wrapper to built in survival functions
//'
//' The function wraps all builtin-survival distributions. User can choose
//' between the \code{weibull}, \code{gamma}, \code{gengamma}(generalized gamma)
//' and \code{burr} (Burr type XII distribution). It is the user responsibility
//' to pass the appropriate list of parameters as follows:
//' \describe{
//' \item{weibull}{\code{scale} (the scale) and \code{shape} (the shape)
//'     parameters.}
//' \item{burr}{\code{scale} (the scale) and \code{shape1} (the shape1) and
//'     \code{shape2} (the shape2) parameters.} 
//' \item{gamma}{ \code{scale} (the scale) and \code{shape} (the shape)
//'     parameter.}
//' \item{gengamma}{\code{mu} (location), \code{sigma} (scale) and \code{Q}
//'     (shape) parameters.}
//' }
//' @param t double, time point where the survival is to be evaluated at.
//' @param distPars \code{Rcpp::List} with distribution specific slots,
//'     see details.
//' @param dist character name of the built-in distribution, see details.
//' @return a double giving the value of the survival function at time point
//' \code{t} at the parameters' values.
//'
//' @examples
//' tt <- 2.5
//' ## weibull
//'
//' distP <- list(scale = 1.2, shape = 1.16)
//' alpha <- exp(-log(distP[["scale"]]) / distP[["shape"]])
//' pweibull(q = tt, scale = alpha, shape = distP[["shape"]],
//'                       lower.tail = FALSE)
//' surv(tt, distP, "weibull") ## (almost) same
//'
//' ## gamma
//' distP <- list(shape = 0.5, rate = 1.0 / 0.7)
//' pgamma(q = tt, rate = distP[["rate"]], shape = distP[["shape"]],
//'                     lower.tail = FALSE)
//' surv(tt, distP, "gamma")  ## (almost) same
//'
//' ## generalized gamma
//' distP <- list(mu = 0.5, sigma = 0.7, Q = 0.7)
//' flexsurv::pgengamma(q = tt, mu = distP[["mu"]],
//'                     sigma = distP[["sigma"]],
//'                     Q = distP[["Q"]],
//'                     lower.tail = FALSE)
//' surv(tt, distP, "gengamma")  ## (almost) same
//'
//' @export
// [[Rcpp::export]]
double surv (double t, const Rcpp::List distPars, const std::string dist) {
  if (dist == "weibull")
    return(sWeibull(t, distPars));
  else if (dist == "gamma")
    return(sgamma(t, distPars));
  else if (dist == "gengamma")
    return(sgengamma(t, distPars));
  else if (dist == "burr")
    return(sBurr(t, distPars));
  else
    stop("distribution not supported !");

  return(0.0); // doesn't each this line but the compiler complains
}

arma::vec getextrapolPars(const Rcpp::List distPars, const std::string dist) {
  arma::vec res(2);

  if (dist == "weibull") {
    res(0) = as <double> (distPars["shape"]) + 1.0;
    res(1) = 2.0;
  } else if (dist == "gamma") {
    res(0) = as <double> (distPars["shape"]) + 1.0;
    res(1) = 2.0;
  } else if (dist == "gengamma") {
    double sig = as <double> (distPars["sigma"]);
    double Q = as <double> (distPars["Q"]);
    res(0) = sig * Q + 1.0;
    res(1) = 2.0;
  } else if (dist == "burr") {
    res(0) = as <double> (distPars["shape1"]) + 1.0;
    res(1) = 2.0;
  } else
    stop("distribution not supported !");

  return(res);
}
