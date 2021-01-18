#' Compute the tau2 estimator of scale
#'
#' @param x A vector of residuals.
#' @return The tau2 estimate of scale.
#' @description The tau2-estimator is a robust measure of the scale. The exact formula of the estimator is in Crevits and Croux (2016), equation 3.10.
#' @references Crevits, R., and Croux, C (2016) "Forecasting with Robust Exponential Smoothing with Damped Trend and Seasonal Components".\emph{Working paper}. \url{https://doi.org/10.13140/RG.2.2.11791.18080}
#' 
#' @examples
#' set.seed(100)
#' e <- 10*rnorm(100)
#' mse <- mean(e^2) 
#' tse <- tau2(e) 
#' @export
tau2 = function(x){
  .Call('robets_tau2', PACKAGE = 'robets', x)
}
