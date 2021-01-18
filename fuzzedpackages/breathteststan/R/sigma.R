#' @title S3 method to exctract the residual standard deviation
#' @description Functions for S3 method defined in breathtestcore for
#' \code{stan_fit} and \code{stan_group fit}.
#' @param object A Stan-based fit
#' @param ... Not used
#' @return A numeric value giving the sigma (= average residual standard deviation) term
#' from the Stan fit.
#' @importFrom stats sigma
#' @export
sigma.breathteststanfit = function(object, ...){
  mean(rstan::extract(object$stan_fit, permuted = TRUE, pars = c( "sigma"))$sigma)
}

