#' @title Simulate Data from Gaussian Graphical Model
#' @description Produce one or more samples from the specified Gaussian graphical model.
#'
#' @param n The number of samples required.
#' @param Omega The inverse covariance matrix of the specified Gaussian graphical model.
#'
#' @return A numeric matrix with \code{n} rows and \code{p} variables where \code{p} corresponds to the dimension of \code{Omega}.
#'
#' @note \code{Omega} should be positive definite.
#'
#' @export
#'
#' @examples
#' library(gif)
#'
#' set.seed(1)
#' n <- 200
#' p <- 100
#' Omega <- diag(1, p, p)
#' for(i in 1:(p - 1)) {
#'   Omega[i, i + 1] <- 0.5
#'   Omega[i + 1, i] <- 0.5
#' }
#' x <- ggm.generator(n, Omega)
ggm.generator <- function(n, Omega) {
  stopifnot(!missing(n) & !missing(Omega))
  stopifnot(n > 0 & round(n, 0) == n)
  stopifnot(isSymmetric(Omega))
  stopifnot(det(Omega) > 0)

  p <- ncol(Omega)
  sigma <- solve(Omega)
  x <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigma)

  return(x)
}
