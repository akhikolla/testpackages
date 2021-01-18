#' Hajek-Rosen variance estimator
#'
#' @description 
#' 
#' Estimator of the variance of the Horvitz-Thompson estimator. It is based on the variance estimator of the conditional Poisson sampling design. See Tillé (2020, Chapter 5) for more informations.
#' 
#' @param y vector of size \eqn{n} that represent the variable of interest.
#' @param pik vector of the inclusion probabilities. The length should be equal to \eqn{n}.
#' @param s index vector of size \eqn{n} with elements equal to the selected units.
#'
#' @details
#' 
#' The function computes the following quantity :
#' 
#' \deqn{v_{HAJ}(\widehat{Y}_{HT}) = \frac{n}{n-1} \sum_{k\in S} (1-\pi_k)\left( \frac{y_k}{\pi_k}-\frac{ \sum_{l\in S} (1-\pi_k)/\pi_k }{\sum_{l\in S} (1-\pi_k) } \right)^2 .}
#'  
#' This estimator is well-defined for maximum entropy sampling design and use only inclusion probabilities of order one.
#' 
#' @return A number, the variance
#' 
#' @references 
#' Tillé, Y. (2020). Sampling and estimation from finite populations. New York: Wiley
#' 
#' @encoding UTF-8
#' @import Matrix
#' 
#' @useDynLib WaveSampling, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#' 
varHAJ <- function(y,pik,s){
  n <- length(s)
  return((n/(n-1))*sum((1-pik)*(y/pik-sum(y*(1-pik)/pik)/sum(1-pik))^2));
}