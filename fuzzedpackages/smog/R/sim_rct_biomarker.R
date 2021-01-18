#' Simulate a randomized clinical trial with biomarkers
#' 
#' \code{sim_rct_biomarker} is used to simulate clinical trial data with
#' specified treatment, prognostic, and predictive effect sizes.
#' 
#' @param n Number of subjects.
#' @param p Number of biomarkers.
#' @param p_prog Number of biomarkers with prognostic effects only.
#' @param p_pred Number of biomarkers with predictive effects only.
#' @param p_both Number of biomarkers with both prognostic and predictive effects
#' @param v_trt Variance of response due to treatment.
#' @param v_prog Variance of response due to prognostic effects.
#' @param v_pred Variance of response due to predictive effects.
#' @param v_err Variance of response due to random noise.
#' @param corr Autocorrelation parameter between biomarkers, default is \code{NULL}.
#' @param family The distribution family for response variable, can be ``gaussian'',
#'               or ``binomial''. Default is ``gaussian''.
#' @param ... further arguments passed to or from other methods. 
#' 
#' @return A list containing several variables.
#'        \describe{
#'        \item{T}{Treatment status in 1 or -1 values.}
#'        \item{X}{Biomarkers.}
#'        \item{W}{Hadamard product of treatment and biomarkers.}
#'        \item{M}{Model matrix - binding of \code{T}, \code{X}, and \code{W}.}
#'        \item{Y}{Response.}
#'        \item{Y0}{Response without error.}
#'        \item{tau}{Treatment effect.}
#'        \item{beta}{Prognostic effects.}
#'        \item{gamma}{Predictive effects.}
#'        \item{theta}{All effects corresponding to \code{M}.}
#'        }
#' @examples
#' sim <- sim_rct_biomarker(n = 1e3)
#' var(as.vector(sim$T * sim$tau))
#' var(as.vector(sim$X %*% sim$beta))
#' var(as.vector(sim$W %*% sim$gamma))
#' 
#' @author Chong Ma \email{chong.ma@@yale.edu}, Kevin Galinsky \email{Kevin.Galinsky@@takeda.com}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @importFrom stats rnorm
#' @export
sim_rct_biomarker <- function (
  n = 50, p = 100,
  p_prog = 5, p_pred = 5, p_both = 5,
  v_trt = 0.4, v_prog = 0.2, v_pred = 0.2, v_err = 0.2,
  corr = NULL, family = "gaussian",...) 
  {
  n_pla <- floor(n / 2)
  n_trt <- n - n_pla
  p1 <- p_prog+p_pred+p_both
  
  # create return list
  sim <- list()

  # set up model matrices
  sim[["T"]] <- c(rep(-1, n_pla), rep(1, n_trt))
  sim[["X"]] <- matrix(rnorm(n * p), nrow = n)
  if(!is.null(corr)){
    S <- chol(corr^abs(row(diag(p1))-col(diag(p1))))
    sim[["X"]][,1:p1] <- sim[["X"]][,1:p1]%*%S
  }
  sim[["W"]] <- sim$T * sim$X
  sim[["M"]] <- cbind(sim$T, sim$X, sim$W)

  # set up effect sizes

  # tau
  # - var(T) = 1
  # - one variable
  sim[["tau"]] <- sqrt(v_trt)

  # - var(X) = 1
  # - p_prog and p_both variables
  beta <- sqrt(v_prog / (p_prog + p_both))

  # - var(W) = E(Var(W|T)) + Var(E(W|T)) = 0.5 [ Var(W|T=-1) + Var(W|T=1) ] = 1
  # - p_pred and p_both variables
  gamma <- sqrt(v_pred / (p_pred + p_both))

  # set up effects
  p_nz <- p_prog + p_pred + p_both
  sim[["beta"]]  <- c(rep(beta, p_prog), rep(0, p_pred), rep(beta, p_both), rep(0, p - p_nz))
  sim[["gamma"]] <- c(rep(0, p_prog), rep(gamma, p_pred), rep(gamma, p_both), rep(0, p - p_nz))

  sim[["theta"]] <- c(sim$tau, sim$beta, sim$gamma)

  sim[["Y0"]] <- sim$M %*% sim$theta
  sim[["Y"]] <- rnorm(n, sim$Y0, sqrt(v_err))
  
  if(family == "binomial"){
    prob = exp(sim[["Y"]])/(1+exp(sim[["Y"]]))
    sim[["Y"]] = ifelse(prob<0.5,0,1)
  }

  return(sim)
}
