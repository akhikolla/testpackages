

#' @title Tangent von Mises--Fisher distribution
#'
#' @description Density and simulation of the Tangent von Mises--Fisher (TM)
#' distribution on
#' \eqn{S^{p-1}:=\{\mathbf{x}\in R^p:||\mathbf{x}||=1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p\ge 2}. The distribution arises
#' by considering the
#' \link[=tang-norm-decomp]{tangent-normal decomposition} with
#' multivariate \link[=cosines-signs]{signs} distributed as a
#' \link[=vMF]{von Mises--Fisher} distribution.
#'
#' @inheritParams tang-norm-decomp
#' @param mu the directional mean \eqn{\boldsymbol{\mu}}{\mu} of the vMF
#' used in the multivariate signs. A unit-norm vector of length \code{p - 1}.
#' @param kappa concentration parameter \eqn{\kappa} of the vMF used in the
#' multivariate signs. A nonnegative scalar.
#' @return
#' Depending on the function:
#' \itemize{
#' \item \code{d_TM}: a vector of length \code{nx} or \code{1} with the
#' evaluated density at \code{x}.
#' \item \code{r_TM}: a matrix of size \code{c(n, p)} with the random sample.
#' }
#' @details
#' The functions are wrappers for \code{\link{d_tang_norm}} and
#' \code{\link{r_tang_norm}} with \code{d_U = \link{d_vMF}} and
#' \code{r_U = \link{r_vMF}}.
#' @author Eduardo García-Portugués, Davy Paindaveine, and Thomas Verdebout.
#' @references
#' García-Portugués, E., Paindaveine, D., Verdebout, T. (2020) On optimal tests
#' for rotational symmetry against new classes of hyperspherical distributions.
#' \emph{Journal of the American Statistical Association}, to appear.
#' \url{https://doi.org/10.1080/01621459.2019.1665527}
#' @examples
#' ## Simulation and density evaluation for p = 2
#'
#' # Parameters
#' p <- 2
#' n <- 1e3
#' theta <- c(rep(0, p - 1), 1)
#' mu <- c(rep(0, p - 2), 1)
#' kappa <- 1
#' kappa_V <- 2
#'
#' # Required functions
#' r_V <- function(n) r_g_vMF(n = n, p = p, kappa = kappa_V)
#' g_scaled <- function(t, log) {
#'   g_vMF(t, p = p - 1, kappa = kappa_V, scaled = TRUE, log = log)
#' }
#'
#' # Sample and color according to density
#' x <- r_TM(n = n, theta = theta, r_V = r_V, mu = 1, kappa = kappa)
#' col <- viridisLite::viridis(n)
#' r <- runif(n, 0.95, 1.05) # Radius perturbation to improve visualization
#' dens <- d_TM(x = x, theta = theta, g_scaled = g_scaled, mu = mu,
#'              kappa = kappa)
#' plot(r * x, pch = 16, col = col[rank(dens)])
#'
#' ## Simulation and density evaluation for p = 3
#'
#' # Parameters
#' p <- 3
#' n <- 5e3
#' theta <- c(rep(0, p - 1), 1)
#' mu <- c(rep(0, p - 2), 1)
#' kappa <- 1
#' kappa_V <- 2
#'
#' # Sample and color according to density
#' x <- r_TM(n = n, theta = theta, r_V = r_V, mu = mu, kappa = kappa)
#' col <- viridisLite::viridis(n)
#' dens <- d_TM(x = x, theta = theta, g_scaled = g_scaled, mu = mu,
#'              kappa = kappa)
#' rgl::plot3d(x, col = col[rank(dens)], size = 5)
#'
#' ## A non-vMF angular function: g(t) = 1 - t^2. It is sssociated to the
#' ## Beta(1/2, (p + 1)/2) distribution.
#'
#' # Scaled angular function
#' g_scaled <- function(t, log) {
#'   log_c_g <- lgamma(0.5 * p) + log(0.5 * p / (p - 1)) - 0.5 * p * log(pi)
#'   log_g <- log_c_g + log(1 - t^2)
#'   switch(log + 1, exp(log_g), log_g)
#' }
#'
#' # Simulation
#' r_V <- function(n) {
#'   sample(x = c(-1, 1), size = n, replace = TRUE) *
#'     sqrt(rbeta(n = n, shape1 = 0.5, shape2 = 0.5 * (p + 1)))
#' }
#'
#' # Sample and color according to density
#' kappa <- 0.5
#' x <- r_TM(n = n, theta = theta, r_V = r_V, mu = mu, kappa = kappa)
#' col <- viridisLite::viridis(n)
#' dens <- d_TM(x = x, theta = theta, g_scaled = g_scaled,
#'              mu = mu, kappa = kappa)
#' rgl::plot3d(x, col = col[rank(dens)], size = 5)
#' @seealso \code{\link{tang-norm-decomp}},
#' \code{\link{tangent-elliptical}}, \code{\link{vMF}}.
#' @name tangent-vMF
#' @aliases TM


#' @rdname tangent-vMF
#' @export
d_TM <- function(x, theta, g_scaled, d_V, mu, kappa, log = FALSE) {

  # Check coherence of theta and mu
  if (length(theta) != (length(mu) + 1)) {

    stop("theta and mu do not have lengths p and p - 1, respectively.")

  }

  # Log-density of the vMF
  d_U <- function(x, log) d_vMF(x = x, mu = mu, kappa = kappa, log = log)

  # Density
  log_dens <- d_tang_norm(x = x, theta = theta, g_scaled = g_scaled,
                          d_V = d_V, d_U = d_U, log = TRUE)
  return(switch(log + 1, exp(log_dens), log_dens))

}


#' @rdname tangent-vMF
#' @export
r_TM <- function(n, theta, r_V, mu, kappa) {

  # Check coherence of theta and mu
  if (length(theta) != (length(mu) + 1)) {

    stop("theta and mu do not have lengths p and p - 1, respectively.")

  }

  # Simulate from a vMF
  r_U <- function(x) r_vMF(n = n, mu = mu, kappa = kappa)

  # Sample
  return(r_tang_norm(n = n, theta = theta, r_V = r_V, r_U = r_U))

}
