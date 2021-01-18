

#' @title Distributions based on the tangent-normal decomposition
#'
#' @description Density and simulation of a distribution on
#' \eqn{S^{p-1}:=\{\mathbf{x}\in R^p:||\mathbf{x}||=1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p\ge 2}, obtained by the
#' tangent-normal decomposition. The \emph{tangent-normal decomposition} of
#' the random vector \eqn{\mathbf{X}\in S^{p-1}}{X \in S^{p-1}} is
#' \deqn{V\boldsymbol{\theta} +
#' \sqrt{1 - V^2}\boldsymbol{\Gamma}_{\boldsymbol{\theta}}\mathbf{U}}{
#' V \theta + \sqrt(1 - V^2) \Gamma_\theta U}
#' where \eqn{V := \mathbf{X}'\boldsymbol{\theta}}{V := X'\theta} is a
#' random variable in \eqn{[-1, 1]} (the \emph{cosines} of
#' \eqn{\mathbf{X}}{X}) and
#' \eqn{\mathbf{U} := \boldsymbol{\Gamma}_{\boldsymbol{\theta}}\mathbf{X}/
#' ||\boldsymbol{\Gamma}_{\boldsymbol{\theta}}\mathbf{X}||}{
#' U := \Gamma_\theta X / ||\Gamma_\theta X||} is a random vector in
#' \eqn{S^{p-2}} (the \emph{multivariate signs} of \eqn{\mathbf{X}}{X})
#' and \eqn{\boldsymbol{\Gamma}_{\boldsymbol{\theta}}}{\Gamma_\theta} is the
#' \eqn{p\times(p-1)}{p x (p-1)} matrix computed by \code{\link{Gamma_theta}}.
#'
#' The tangent-normal decomposition can be employed for constructing
#' distributions for \eqn{\mathbf{X}}{X} that arise for certain choices of
#' \eqn{V} and \eqn{\mathbf{U}}{U}. If \eqn{V} and
#' \eqn{\mathbf{U}}{U} are \emph{independent}, then simulation from
#' \eqn{\mathbf{X}}{X} is straightforward using the tangent-normal
#' decomposition. Also, the density of \eqn{\mathbf{X}}{X} at
#' \eqn{\mathbf{x}\in S^{p-1}}{x \in S^{p-1}},
#' \eqn{f_\mathbf{X}(\mathbf{x})}{f_X(x)}, is readily computed as
#' \deqn{f_\mathbf{X}(\mathbf{x})=
#' \omega_{p-1}c_g g(t)(1-t^2)^{(p-3)/2}f_\mathbf{U}(\mathbf{u})}{
#' f_X(x) = \omega_{p-1} c_g g(t) (1-t^2)^{(p-3)/2} f_U(u)}
#' where \eqn{t:=\mathbf{x}'\boldsymbol{\theta}}{v := x'\theta},
#' \eqn{\mathbf{u}:=\boldsymbol{\Gamma}_{\boldsymbol{\theta}}\mathbf{x}/
#' ||\boldsymbol{\Gamma}_{\boldsymbol{\theta}}\mathbf{x}||}{
#' u := \Gamma_\theta x / ||\Gamma_\theta x||},
#' \eqn{f_\mathbf{U}}{f_U} is the density of \eqn{\mathbf{U}}{U},
#' and \eqn{f_V(v) := \omega_{p-1} c_g g(v) (1 - v^2)^{(p-3)/2}} is the density
#' of \eqn{V} for an angular function \eqn{g} with normalizing constant
#' \eqn{c_g}. \eqn{\omega_{p-1}} is the surface area of \eqn{S^{p-2}}.
#'
#' @inheritParams unif
#' @param theta a unit norm vector of size \code{p} giving the axis of
#' rotational symmetry.
#' @param g_scaled the \emph{scaled} angular density \eqn{c_g g}. In the
#' form \cr\code{g_scaled <- function(t, log = TRUE) {...}}. See examples.
#' @param d_V the density \eqn{f_V}. In the form
#' \code{d_V <- function(v, log = TRUE) {...}}. See examples.
#' @param d_U the density \eqn{f_\mathbf{U}}{f_U}. In the form
#' \code{d_U <- function(u, log = TRUE) {...}}. See examples.
#' @param r_U a function for simulating \eqn{\mathbf{U}}{U}. Its first argument
#' must be the sample size. See examples.
#' @param r_V a function for simulating \eqn{V}. Its first argument must be
#' the sample size. See examples.
#' @return
#' Depending on the function:
#' \itemize{
#' \item \code{d_tang_norm}: a vector of length \code{nx} or \code{1}
#' with the evaluated density at \code{x}.
#' \item \code{r_tang_norm}: a matrix of size \code{c(n, p)} with the
#' random sample.
#' }
#' @details
#' Either \code{g_scaled} or \code{d_V} can be supplied to \code{d_tang_norm}
#' (the rest of the arguments are compulsory). One possible choice for
#' \code{g_scaled} is \code{\link{g_vMF}} with \code{scaled = TRUE}. Another
#' possible choice is the angular function \eqn{g(t) = 1 - t^2}, normalized by
#' its normalizing constant
#' \eqn{c_g = (\Gamma(p/2) p) / (2\pi^{p/2} (p - 1))} (see examples).
#' This angular function makes \eqn{V^2} to be distributed as a
#' \eqn{\mathrm{Beta}(1/2,(p+1)/2)}{Beta(1/2, (p+1)/2)}.
#'
#' The normalizing constants and densities are computed through log-scales for
#' numerical accuracy.
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
#' n <- 1e3
#' p <- 2
#' theta <- c(rep(0, p - 1), 1)
#' mu <- c(rep(0, p - 2), 1)
#' kappa_V <- 2
#' kappa_U <- 0.1
#'
#' # The vMF scaled angular function
#' g_scaled <- function(t, log) {
#'   g_vMF(t, p = p - 1, kappa = kappa_V, scaled = TRUE, log = log)
#' }
#'
#' # Cosine density for the vMF distribution
#' d_V <- function(v, log) {
#'  log_dens <- g_scaled(v, log = log) + (p - 3)/2 * log(1 - v^2)
#'  switch(log + 1, exp(log_dens), log_dens)
#' }
#'
#' # Multivariate signs density based on a vMF
#' d_U <- function(x, log) d_vMF(x = x, mu = mu, kappa = kappa_U, log = log)
#'
#' # Simulation functions
#' r_V <- function(n) r_g_vMF(n = n, p = p, kappa = kappa_V)
#' r_U <- function(n) r_vMF(n = n, mu = mu, kappa = kappa_U)
#'
#' # Sample and color according to density
#' x <- r_tang_norm(n = n, theta = theta, r_V = r_V, r_U = r_U)
#' r <- runif(n, 0.95, 1.05) # Radius perturbation to improve visualization
#' col <- viridisLite::viridis(n)
#' dens <- d_tang_norm(x = x, theta = theta, g_scaled = g_scaled, d_U = d_U)
#' # dens <- d_tang_norm(x = x, theta = theta, d_V = d_V, d_U = d_U) # The same
#' plot(r * x, pch = 16, col = col[rank(dens)])
#'
#' ## Simulation and density evaluation for p = 3
#'
#' # Parameters
#' p <- 3
#' n <- 5e3
#' theta <- c(rep(0, p - 1), 1)
#' mu <- c(rep(0, p - 2), 1)
#' kappa_V <- 2
#' kappa_U <- 2
#'
#' # Sample and color according to density
#' x <- r_tang_norm(n = n, theta = theta, r_V = r_V, r_U = r_U)
#' col <- viridisLite::viridis(n)
#' dens <- d_tang_norm(x = x, theta = theta, g_scaled = g_scaled, d_U = d_U)
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
#' # Cosine density
#' d_V <- function(v, log) {
#'   log_dens <- w_p(p = p - 1, log = TRUE) + g_scaled(t = v, log = TRUE) +
#'     (0.5 * (p - 3)) * log(1 - v^2)
#'   switch(log + 1, exp(log_dens), log_dens)
#' }
#'
#' # Simulation
#' r_V <- function(n) {
#'   sample(x = c(-1, 1), size = n, replace = TRUE) *
#'     sqrt(rbeta(n = n, shape1 = 0.5, shape2 = 0.5 * (p + 1)))
#' }
#'
#' # Sample and color according to density
#' r_U <- function(n) r_unif_sphere(n = n, p = p - 1)
#' x <- r_tang_norm(n = n, theta = theta, r_V = r_V, r_U = r_U)
#' col <- viridisLite::viridis(n)
#' dens <- d_tang_norm(x = x, theta = theta, d_V = d_V, d_U = d_unif_sphere)
#' # dens <- d_tang_norm(x = x, theta = theta, g_scaled = g_scaled,
#' #                     d_U = d_unif_sphere) # The same
#' rgl::plot3d(x, col = col[rank(dens)], size = 5)
#' @seealso \code{\link{Gamma_theta}}, \code{\link{signs}},
#' \code{\link{tangent-elliptical}}, \code{\link{tangent-vMF}},
#' \code{\link{vMF}}.
#' @name tang-norm-decomp


#' @rdname tang-norm-decomp
#' @export
d_tang_norm <- function(x, theta, g_scaled, d_V, d_U, log = FALSE) {

  # x as matrix
  if (is.null(dim(x))) {

    x <- rbind(x)

  }

  # Detect edge cases
  x <- check_unit_norm(x = x, warnings = TRUE)
  theta <- check_unit_norm(x = theta, warnings = TRUE)

  # Dimension
  p <- ncol(x)
  if (p < 2) {

    stop("p must be an integer larger or equal to 2.")

  }
  if (p != length(theta)) {

    stop("x and theta do not have the same dimension.")

  }

  # Compute signs and cosines
  u_theta_x <- signs(X = x, theta = theta, check_X = FALSE)
  v_theta_x <- cosines(X = x, theta = theta, check_X = FALSE)

  # Compute density of V
  if (!missing(d_V)) {

    log_d_V_theta_x <- d_V(v_theta_x, log = TRUE) -
      0.5 * (p - 3) * log(1 - v_theta_x * v_theta_x)

  } else if (!missing(g_scaled)) {

    log_d_V_theta_x <- w_p(p = p - 1, log = TRUE) +
      g_scaled(v_theta_x, log = TRUE)

  } else {

    stop("Must provide g_scaled or d_V")

  }

  # Density by the tangent-normal decomposition
  log_dens <- log_d_V_theta_x + d_U(u_theta_x, log = TRUE)
  return(switch(log + 1, exp(log_dens), log_dens))

}


#' @rdname tang-norm-decomp
#' @export
r_tang_norm <- function(n, theta, r_U, r_V) {

  # Detect edge cases
  theta <- check_unit_norm(x = theta, warnings = TRUE)

  # Only if valid for p >= 2
  p <- length(theta)
  if (p < 2) {

    stop("p must be an integer larger or equal to 2.")

  }

  # Sample cosines
  V <- drop(r_V(n))

  # Sample multivariate signs
  U <- r_U(n) %*% t(Gamma_theta(theta = theta))

  # Sample by the tangent-normal decomposition
  return(V * matrix(theta, nrow = n, ncol = p, byrow = TRUE) +
           sqrt(1 - V * V) * U)

}
