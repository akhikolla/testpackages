

#---------------------------------------------
ypbp.mle <- function(status, Z, degree, tau, g, G,
                     baseline=c("hazard", "odds"), hessian, ...) {

  n <- nrow(Z)
  q <- ncol(Z)
  baseline <- match.arg(baseline)
  if(baseline=="hazard"){
    M <- 1
  }else{
    M <- 2
  }

  hyper_parms = list(h1_gamma=0, h2_gamma=4,
                     mu_psi=0, sigma_psi=4,
                     mu_phi=0, sigma_phi=4)

  stan_data <- list(status=status, Z=Z, q=q, n=n,
                    m=degree, M=M, approach=0, tau=tau,
                    g=g, G=G,
                    h1_gamma=hyper_parms$h1_gamma,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$mu_phi,
                    h2_gamma=hyper_parms$h2_gamma,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_phi)

  fit <- rstan::optimizing(stanmodels$ypbp, data=stan_data,
                           hessian=hessian, ...)

  fit$par <- fit$par[-grep("loglik", names(fit$par))]
  #fit$par <- fit$par[-grep("log_gamma", names(fit$par))]
  fit$theta_tilde <- fit$theta_tilde[-grep("loglik", names(fit$theta_tilde))]

  return(fit)
}

#---------------------------------------------
ypbp.bayes <- function(status, Z, degree, tau, g, G,
                       baseline=c("hazard", "odds"), hyper_parms, ...) {

  n <- nrow(Z)
  q <- ncol(Z)
  baseline <- match.arg(baseline)
  if(baseline=="hazard"){
    M <- 1
  }else{
    M <- 2
  }

  stan_data <- list(status=status, Z=Z, q=q, n=n,
                    m=degree, M=M, approach=1, tau=tau,
                    g=g, G=G,
                    h1_gamma=hyper_parms$h1_gamma,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$mu_phi,
                    h2_gamma=hyper_parms$h2_gamma,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_phi)

  pars <- c("psi", "phi", "gamma", "loglik")
  fit <- rstan::sampling(stanmodels$ypbp, data=stan_data, pars=pars, ...)

  return(fit)
}

#---------------------------------------------
ypbp2.mle <- function(status, Z, X, degree, tau, g, G,
                      baseline=c("hazard", "odds"), hessian, ...) {

  n <- nrow(Z)
  q <- ncol(Z)
  p <- ncol(X)
  baseline <- match.arg(baseline)
  if(baseline=="hazard"){
    M <- 3
  }else{
    M <- 4
  }

  hyper_parms = list(h1_gamma=0, h2_gamma=4,
                     mu_psi=0, sigma_psi=4,
                     mu_phi=0, sigma_phi=4,
                     mu_beta=0, sigma_beta=4)

  stan_data <- list(status=status, Z=Z, X=X, q=q, p=p,
                    n=n, g=g, G=G,
                    m=degree, M=M, approach=0, tau=tau,
                    h1_gamma=hyper_parms$h1_gamma,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$mu_phi,
                    mu_beta=hyper_parms$mu_beta,
                    h2_gamma=hyper_parms$h2_gamma,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_phi,
                    sigma_beta=hyper_parms$sigma_beta)

  fit <- rstan::optimizing(stanmodels$ypbp2, data=stan_data,
                           hessian=hessian, ...)

  fit$par <- fit$par[-grep("loglik", names(fit$par))]
  fit$theta_tilde <- fit$theta_tilde[-grep("loglik", names(fit$theta_tilde))]

  return(fit)
}

#---------------------------------------------
ypbp2.bayes <- function(status, Z, X, degree, tau,
                        g, G, baseline=c("hazard", "odds"),
                        hyper_parms, ...) {

  n <- nrow(Z)
  q <- ncol(Z)
  p <- ncol(X)
  baseline <- match.arg(baseline)
  if(baseline=="hazard"){
    M <- 3
  }else{
    M <- 4
  }

  stan_data <- list(status=status, Z=Z, X=X, q=q, p=p, n=n,
                    m=degree, M=M, approach=1, tau=tau, g=g, G=G,
                    h1_gamma=hyper_parms$h1_gamma,
                    mu_psi=hyper_parms$mu_psi,
                    mu_phi=hyper_parms$mu_phi,
                    mu_beta=hyper_parms$mu_beta,
                    h2_gamma=hyper_parms$h2_gamma,
                    sigma_psi=hyper_parms$sigma_psi,
                    sigma_phi=hyper_parms$sigma_phi,
                    sigma_beta=hyper_parms$sigma_beta)

  fit <- rstan::sampling(stanmodels$ypbp2, data=stan_data, ...)

  return(fit)
}

#---------------------------------------------

#' Fits the Yang and Prentice using Bernstein polynomials to model the baseline distribution.
#' @aliases{ypbp}
#' @export
#' @description Fits the Yang and Prentice model with either the baseline hazard hazard or the baseline odds modeled via Bernstein polynomials.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which ypbp is called.
#' @param degree number of intervals of the PE distribution. If NULL, default value (square root of n) is used.
#' @param tau the maximum time of follow-up. If NULL, tau = max(time), where time is the vector of observed survival times.
#' @param approach approach to be used to fit the model (mle: maximum likelihood; bayes: Bayesian approach).
#' @param baseline baseline function to be modeled.
#' @param hessian logical; If TRUE (default), the hessian matrix is returned when approach="mle".
#' @param hyper_parms a list containing the hyper-parameters of the prior distributions (when approach = "bayes"). If not specified, default values are used.
#' @param ... Arguments passed to either `rstan::optimizing` or `rstan::sampling` .
#' @return ypbp returns an object of class "ypbp" containing the fitted model.
#'
#' @examples
#' \donttest{
#' library(YPBP)
#' mle1 <- ypbp(Surv(time, status)~trt, data=gastric, baseline = "hazard")
#' mle2 <- ypbp(Surv(time, status)~trt, data=gastric, baseline = "odds")
#' bayes1 <- ypbp(Surv(time, status)~trt, data=gastric, baseline = "hazard",
#'                approach = "bayes", chains = 2, iter = 500)
#' bayes2 <- ypbp(Surv(time, status)~trt, data=gastric, baseline = "odds",
#'                approach = "bayes", chains = 2, iter = 500)
#' }
#'
#'
ypbp <- function(formula, data, degree=NULL, tau=NULL,
                 approach = c("mle", "bayes"), baseline=c("hazard", "odds"),
                 hessian=TRUE, hyper_parms = list(h1_gamma=0, h2_gamma=4,
                                    mu_psi=0, sigma_psi=4,
                                    mu_phi=0, sigma_phi=4,
                                    mu_beta=0, sigma_beta=4), ...){

  approach <- match.arg(approach)
  baseline <- match.arg(baseline)
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  resp <- stats::model.response(mf)
  time <- resp[,1]
  status <- resp[,2]
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- suppressWarnings(try( stats::model.matrix(formula, data = mf, rhs = 2), TRUE))
  labels <- colnames(Z)[-1]
  labels.ph <- colnames(X)[-1]
  Z <- matrix(Z[,-1], ncol=length(labels))
  if(ncol(X)>0){
    labels.ph <- colnames(X)[-1]
    X <- matrix(X[,-1], ncol=length(labels.ph))
  }

  n <- nrow(Z)
  q <- ncol(Z)
  p <- ncol(X)


  if(is.null(tau)){
    tau <- max(time)
  }

  if(is.null(degree)){
    degree <- ceiling(sqrt(length(time)))
  }

  bases <- bp(time, degree, tau)
  g <- bases$b
  G <- bases$B

  if(approach=="mle"){
    if(p==0){
      fit <- ypbp.mle(status=status, Z=Z,
                      degree=degree, tau=tau, g=g, G=G,
                      baseline=baseline, hessian=hessian, ...)
    }else{
      fit <- ypbp2.mle(status=status, Z=Z, X=X,
                       degree=degree, tau=tau, g=g, G=G,
                       baseline=baseline, hessian=hessian, ...)
    }
  }else{
    if(p==0){
      fit <- ypbp.bayes(status=status, Z=Z,
                        degree=degree, tau=tau, g=g, G=G,
                        baseline=baseline, hyper_parms=hyper_parms, ...)
    }else{
      fit <- ypbp2.bayes(status=status, Z=Z, X=X,
                         degree=degree, tau=tau, g=g, G=G,
                         baseline=baseline, hyper_parms=hyper_parms, ...)
    }
  }

  output <- list(fit=fit)

  output$n <- n
  output$q <- q
  output$p <- p

  output$degree <- degree
  output$tau <- tau
  output$call <- match.call()
  output$formula <- formula
  output$terms <- stats::terms.formula(formula)
  output$mf <- mf
  output$labels <- labels
  output$approach <- approach
  output$baseline <- baseline

  if(p>0){
    output$labels.ph <- labels.ph
  }


  class(output) <- "ypbp"
  return(output)
}


#---------------------------------------------

ypbpBoot <- function(formula, data, degree=NULL, tau=NULL,
                     nboot = 4000, ...){
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  resp <- stats::model.response(mf)
  time <- resp[,1]
  status <- resp[,2]
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- suppressWarnings(try( stats::model.matrix(formula, data = mf, rhs = 2), TRUE))
  labels <- colnames(Z)[-1]
  labels.ph <- colnames(X)[-1]
  Z <- matrix(Z[,-1], ncol=length(labels))
  if(ncol(X)>0){
    labels.ph <- colnames(X)[-1]
    X <- matrix(X[,-1], ncol=length(labels.ph))
  }

  n <- nrow(Z)
  q <- ncol(Z)
  p <- ncol(X)


  if(is.null(tau)){
    tau <- max(time)
  }

  if(is.null(degree)){
    degree <- ceiling(sqrt(length(time)))
  }


  index <- 1:n
  index1 <- which(status==1)
  index2 <- which(status==0)
  n1 <- length(index1)
  n2 <- length(index2)
  par <- matrix(nrow=nboot, ncol=(2*q+p+degree))

  for(step in 1:nboot){
    samp1 <- sample(index1, size=n1, replace=TRUE)
    samp2 <- sample(index2, size=n2, replace=TRUE)
    samp <- c(samp1, samp2)
    suppressWarnings({invisible(utils::capture.output(object <- ypbp(formula, data=data[samp,], degree=degree, tau=tau, hessian=FALSE, approach="mle", init=0)))})
    if(class(object)!="try-error"){
      par[step, ] <- object$fit$par
      step <- step + 1
    }
  }

  colnames(par) <- names(object$fit$par[-grep("log_", names(object$fit$par))])

  return(par)
}
