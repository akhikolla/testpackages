#' Flexible Regression Models for Proportions
#'
#' @description The function fits some flexible regression models for proportions via a Bayesian approach to inference based on Hamiltonian Monte Carlo algorithm.
#' Available regression models are the flexible beta regression model (\code{type="FB"}, default), the variance inflated beta (\code{type="VIB"}), and the beta one (\code{type="Beta"}).
#'
#'
#' @param formula an object of class \code{`formula`}: a symbolic description of the model to be fitted (of type \code{y ~ x} or \code{y ~ x | z}).
#' @param data an optional data frame, list, or object coercible to a data frame through \code{base::as.data.frame} containing the variables in the model. If not found in data, the variables in formula are taken from the environment from which the function flexreg is called.
#' @param type a character specifying the type of regression model. Current options are the flexible beta regression model \code{"FB"} (default), the variance inflated beta \code{"VIB"}, and the beta one \code{"Beta"}.
#' @param link.mu a character specifying the link function for the mean model (mu). Currently, \code{"logit"} (default), \code{"probit"}, \code{"cloglog"}, and \code{"loglog"} are supported.
#' @param prior.beta a character specifying the prior distribution for the \code{beta} regression coefficients of the mean model. Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.beta a numeric (vector of length 1) specifying the hyperprior parameter for the prior distribution of \code{beta} regression coefficients. A value of 100 is suggested if the prior is \code{"normal"}, 2.5 if \code{"cauchy"}.
#' @param link.phi a character specifying the link function for the precision model (phi). Currently, \code{"identity"} (default), \code{"log"} and \code{"sqrt"} are supported.
#' @param prior.phi a character specifying the prior distribution for precision parameter \code{phi} if \code{link.phi = "identity"}. Currently, \code{"gamma"} (default) and \code{"unif"} are supported.
#' @param hyperparam.phi a numeric (vector of length 1) specifying the hyperprior parameter for the prior distribution of \code{phi}. A value of 0.001 is suggested if the prior is \code{"gamma"}. If the prior is \code{"uniform"} the hyperparameter must be specified to define the upper limit of the support of \code{phi}.
#' @param prior.psi a character specifying the prior distribution for \code{psi} regression coefficients of the precision model (not supported if link.phi is \code{"identity"}). Currently, \code{"normal"} (default) and \code{"cauchy"} are supported.
#' @param hyperparam.psi a numeric (vector of length 1) specifying the hyperprior parameter for the prior distribution of \code{psi} regression coefficients. A value of 100 is suggested if the prior is \code{"normal"}, 2.5 if \code{"cauchy"}.
#' @param n.iter 	a positive integer specifying the number of iterations for each chain (including warmup). The default is 5000.
#' @param burnin.perc the percentage of iterations (out of \code{n.iter}) per chain to discard.
#' @param n.chain a positive integer specifying the number of Markov chains. The default is 1.
#' @param thin a positive integer specifying the period for saving samples. The default is 1.
#' @param verbose \code{TRUE} (default) or \code{FALSE}: flag indicating whether to print intermediate output.
#' @param ... additional arguments for \code{rstan::sampling}.
#'
#' @return The \code{flexreg} function returns an object of class \code{`flexreg`}, i.e. a list with the following elements:
#' \item{\code{call}}{the function call.}
#' \item{\code{formula}}{the original formula.}
#' \item{\code{link.mu}}{a character specifing the link function in the mean model.}
#' \item{\code{link.phi}}{a character specifing the link function in the precision model.}
#' \item{\code{model}}{an object of class \code{`stanfit`} containing the fitted model.}
#' \item{\code{response}}{the response variable, assuming values in (0, 1).}
#' \item{\code{design.X}}{the design matrix for the mean model.}
#' \item{\code{design.Z}}{the design matrix for the precision model (if defined).}
#'
#' @details Let \eqn{\mu} be the mean of a random variable Y whose distribution can be specified in the \code{type} argument.
#' Then the \code{flexreg} function links the parameter \eqn{\mu} to a linear predictor through a function  \eqn{g(\cdot)} specified in \code{link.mu}:
#' \deqn{g(\mu_i) = x_i^t \bold{\beta},} where \eqn{\bold{\beta}} is the vector of regression coefficients for the mean model.
#' By default, the precision parameter \eqn{\phi} is assumed to be constant.
#' It is possible to extend the model by linking \eqn{\phi} to an additional (possibly overlapping) set of covariates through a proper link
#' function \eqn{q(\cdot)}  specified in the \code{link.phi} argument: \deqn{q(\phi_i) = z_i^t \bold{\psi},} where \eqn{\bold{\psi}} is the vector of regression coefficients for the precision model.
#' In \code{flexreg}, the regression model for the mean and, where appropriate, for the precision parameter can be specified in the
#' \code{formula} argument with a formula of type \eqn{y \sim x_1 + x_2 | z_1 + z_2} where covariates on the left of ("|") are included in the regression model
#' for the mean and covariates on the right of ("|") are included in the regression model for the precision.
#'
#' If the second part is omitted, i.e., \eqn{y \sim x_1 + x_2}, the precision is assumed constant for each observation.
#'
#' @references {
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018) A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079 \cr
#' \cr
#' Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020) Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309. doi:10.1177/1471082X18821213 \cr
#' \cr
#'  Ferrari, S.L.P., and Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. Journal of Applied Statistics, \bold{31}(7), 799--815. doi:10.1080/0266476042000214501
#' }
#'
#' @examples
#' data("Reading")
#' FB <- flexreg(accuracy ~ iq, Reading, type="FB", n.iter=1000)
#'
#' @import Rcpp methods rstan
#'
#' @export

flexreg <- function(formula, data, type="FB", link.mu="logit",
                    prior.beta = "normal", hyperparam.beta = 100,
                    link.phi = NULL, prior.phi = NULL, hyperparam.phi = NULL,
                    prior.psi = NULL, hyperparam.psi = NULL,
                    n.iter=5000, burnin.perc=.5, n.chain=1, thin=1, verbose=TRUE, ...)
  {
  cl <- match.call()

  if(is.na(match(type, c("FB", "Beta", "VIB")))) stop("Please specify the type of model correctly.")

  if (missing(data)) data <- environment(formula) else  data <- as.data.frame(data)#

  formula <- Formula::as.Formula(formula)
  if (is.character(link.phi)) link_code_phi <- pmatch(link.phi, c("identity", "log", "sqrt"))

  if(length(formula)[2] >= 2){
    model.phi <- TRUE
    if(is.null(link.phi)) link_code_phi <- 2
    if(link_code_phi <2) stop("Error: invalid link function for regression model for phi")
    if(!is.null(prior.phi) | !is.null(hyperparam.phi)) stop("Error: prior chosen for phi but regression model for phi in formula. Please define priors for coefficients psi instead.")
    if(length(formula)[2] > 2)  warning("formula must not have more than two RHS parts")
  } else if(length(formula)[2] < 2){
  model.phi <- FALSE
  if(is.null(link.phi)) link_code_phi <- 1
  if(link_code_phi >1) stop("Error: covariates for regression model for phi must be specified")
  if(!is.null(prior.psi) | !is.null(hyperparam.psi)) stop("Error: prior for psi coefficients defined but no regression model for phi. Please define covariates for regression model for phi")
  }

  mf <- model.frame(formula, data)
  y <- model.response(mf, "numeric")
  N <- length(y)

  if(N < 1)
    stop("Empty model")
  if(!(min(y) > 0 & max(y) < 1))
    stop("Invalid dependent variable, all observations must be in (0, 1).")#funzione normalize?

  X <- model.matrix(formula, data = mf, rhs = 1)
  if(isFALSE(model.phi)){
    Z <- NULL} else {
       Z <- model.matrix(formula, data = mf, rhs = 2)
    }

  if(is.character(link.mu)) link_code_mu <- pmatch(link.mu, c("logit", "probit", "cloglog", "loglog")) else
    stop("Invalid link function for mu.")

  if(is.null(prior.beta)) link_prior_beta <- 1 else
    if (is.character(prior.beta)) link_prior_beta  <- pmatch(prior.beta, c("normal", "cauchy")) else
      stop("Invalid prior for beta parameter.")

  if(hyperparam.beta < 0) stop("Hyperprior for beta coefficients must be positive")

  if(link_prior_beta == 1 & hyperparam.beta < 10) warning("An hyperprior of at least 10 for normal prior for beta coefficients is recommended")

  if(link_prior_beta == 2 & hyperparam.beta != 2.5) warning("An hyperprior of 2.5 for cauchy prior for beta coefficients is recommended")
  if(link_code_phi == 1) {
    prior_code_phi <- ifelse(is.null(prior.phi),1,
                             pmatch(prior.phi, c( "gamma", "unif")))
    if(prior_code_phi==1) hyper_phi <-ifelse(is.null(hyperparam.phi),0.001, hyperparam.phi) else #add warnings if g grande
      if(prior_code_phi==2) hyper_phi <-ifelse(is.null(hyperparam.phi), stop("Please specify an hyperparameter A>0 for uniform prior for phi"), hyperparam.phi)
      link_prior_psi <- NULL
  }  else if (link_code_phi != 1) {
    prior_code_phi <- NULL
    hyper_phi <- NULL
    if(is.null(prior.psi)) link_prior_psi <- 1 else
      if (is.character(prior.psi)) link_prior_psi  <- pmatch(prior.psi, c("normal", "cauchy")) else
        stop("Invalid prior for psi parameter.")

    if(is.null(hyperparam.psi)) hyperprior_psi <- ifelse(link_prior_psi==1, 100, 2.5) else
      if(hyperparam.psi < 0) stop("Hyperprior for psi coefficients must be positive") else {
        if(link_prior_psi == 1 & hyperparam.psi < 10) warning("An hyperprior of at least 10 for normal for psi coefficients is recommended")
        if(link_prior_psi == 2 & hyperparam.psi != 2.5) warning("An hyperprior of 2.5 for cauchy prior for psi coefficients is recommended")
      }
    }

    model <- fit.model(model.phi = model.phi, type = type, N = N,  y = y,  X = X, Z = Z, link_code_mu = link_code_mu,
                    link_prior_beta, hyperparam.beta,
                    link_code_phi = link_code_phi,
                    prior_code_phi = prior_code_phi, hyper_phi = hyper_phi,
                    link_prior_psi, hyperparam.psi,
                    n.iter, burnin.perc, n.chain, thin, verbose, ...)

  #questa riga serve per summary
  link.phi <- c("identity", "log", "sqrt")[link_code_phi]
  output <- list(call=cl, formula=formula, link.mu=link.mu, link.phi=link.phi,
                 model=model, response=y, design.X=X, design.Z=Z)
  class(output)<-"flexreg"
  invisible(output)
  #return(output)
}


fit.model <- function(model.phi = model.phi, type = type, N = N,  y = y, X = X, Z = Z, link_code_mu = link_code_mu,
                   link_prior_beta, hyperparam.beta,
                   link_code_phi = link_code_phi,
                   prior_code_phi = prior_code_phi, hyper_phi = hyper_phi,
                   link_prior_psi, hyperparam.psi,
                   n.iter, burnin.perc, n.chain, thin, verbose, ...){

  data.stan <- list(
    N = N,
    y = y,
    X = X,
    Z = Z,
    K = ncol(X),
    H = ncol(Z),
    link_code_mu = link_code_mu,
    link_prior_beta = link_prior_beta,
    hyperprior_beta = hyperparam.beta,
    link_code_phi = link_code_phi,
    prior_code_phi = prior_code_phi,
    hyper_phi = hyper_phi,
    link_prior_psi = link_prior_psi,
    hyperprior_psi = hyperparam.psi
  )

  if(model.phi) stan.model <- stanmodels[[paste(type, "_phi", sep="")]]
  else stan.model <- stanmodels[[type]]

  fit.stan = rstan::sampling(
    object = stan.model,
    data = data.stan,
    chains = n.chain,
    thin = thin,
    #init = initl,
    iter = n.iter, warmup = round(n.iter*burnin.perc),
    refresh = verbose*n.iter/10, #show an update @ each %10
    # seed=1
    ...
  )

  output <- fit.stan
  return(output)

}
