

#' @title von Mises--Fisher distribution
#'
#' @description Density and simulation of the von Mises--Fisher (vMF)
#' distribution on
#' \eqn{S^{p-1}:=\{\mathbf{x}\in R^p:||\mathbf{x}||=1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p\ge 1}. The density at
#' \eqn{\mathbf{x} \in S^{p-1}}{x \in S^{p-1}} is given by
#' \deqn{c^{\mathrm{vMF}}_{p,\kappa}
#' e^{\kappa\mathbf{x}' \boldsymbol{\mu}}
#' \quad\mathrm{with}\quad c^{\mathrm{vMF}}_{p,\kappa}:=
#' \kappa^{(p-2)/2}/((2\pi)^{p/2} I_{(p-2)/2}(\kappa))}{
#' c^{vMF}_{p,\kappa} e^{\kappa x' \mu} with
#' c^{vMF}_{p,\kappa} :=
#' \kappa^{(p-2)/2}/((2\pi)^{p/2} I_{(p-2)/2}(\kappa))}
#' where \eqn{\boldsymbol{\mu}\in S^{p-1}}{\mu \in S^{p-1}} is the directional
#' mean, \eqn{\kappa\ge 0} is the concentration parameter about
#' \eqn{\boldsymbol{\mu}}{\mu}, and \eqn{I_\nu} is the order-\eqn{\nu}
#' modified Bessel function of the first kind.
#'
#' The angular function of the vMF is \eqn{g(t) := e^{\kappa t}}. The
#' associated \emph{cosines} density is
#' \eqn{\tilde g(v):= \omega_{p-1} c^{\mathrm{vMF}}_{p,\kappa}
#' g(v) (1 - v^2)^{(p-3)/2}}{
#' \tilde g(v):= \omega_{p-1} c^{vMF}_{p,\kappa} g(v)(1 - v^2)^{(p-3)/2}},
#' where \eqn{\omega_{p-1}} is the surface area of \eqn{S^{p-2}}.
#'
#' @inheritParams unif
#' @param mu the directional mean \eqn{\boldsymbol{\mu}}{\mu} of the vMF.
#' A unit-norm vector of length \code{p}.
#' @param kappa concentration parameter \eqn{\kappa} of the vMF.
#' A nonnegative scalar. Can be a vector for \code{c_vMF}.
#' @param t a vector with values in \eqn{[-1, 1]}.
#' @param scaled whether to scale the angular function by the von Mises--Fisher
#' normalizing constant. Defaults to \code{TRUE}.
#' @return
#' Depending on the function:
#' \itemize{
#' \item \code{d_vMF}: a vector of length \code{nx} or \code{1} with the
#' evaluated density at \code{x}.
#' \item \code{r_vMF}: a matrix of size \code{c(n, p)} with the random sample.
#' \item \code{c_vMF}: the normalizing constant.
#' \item \code{g_vMF}: a vector of size \code{length(t)} with the evaluated
#' angular function.
#' \item \code{r_g_vMF}: a vector of length \code{n} containing simulated
#' values from the cosines density associated to the angular function.
#' }
#' @details
#' \code{r_g_vMF} implements algorithm VM in Wood (1994). \code{c_vMF} is
#' vectorized on \code{p} and \code{kappa}.
#' @author Eduardo García-Portugués, Davy Paindaveine, and Thomas Verdebout.
#' @references
#' Wood, A. T. A. (1994) Simulation of the von Mises Fisher distribution.
#' \emph{Commun. Stat. Simulat.}, 23(1):157--164.
#' \url{https://doi.org/10.1080/03610919408813161}
#' @examples
#' # Simulation and density evaluation for p = 2
#' mu <- c(0, 1)
#' kappa <- 2
#' n <- 1e3
#' x <- r_vMF(n = n, mu = mu, kappa = kappa)
#' col <- viridisLite::viridis(n)
#' r <- runif(n, 0.95, 1.05) # Radius perturbation to improve visualization
#' plot(r * x, pch = 16, col = col[rank(d_vMF(x = x, mu = mu, kappa = kappa))])
#'
#' # Simulation and density evaluation for p = 3
#' mu <- c(0, 0, 1)
#' kappa <- 2
#' x <- r_vMF(n = n, mu = mu, kappa = kappa)
#' rgl::plot3d(x, col = col[rank(d_vMF(x = x, mu = mu, kappa = kappa))],
#'             size = 5)
#'
#' # Cosines density
#' g_tilde <- function(t, p, kappa) {
#'   exp(w_p(p = p - 1, log = TRUE) +
#'         g_vMF(t = t, p = p, kappa = kappa, scaled = TRUE, log = TRUE) +
#'         ((p - 3) / 2) * log(1 - t^2))
#' }
#'
#' # Simulated data from the cosines density
#' n <- 1e3
#' p <- 3
#' kappa <- 2
#' hist(r_g_vMF(n = n, p = p, kappa = kappa), breaks = seq(-1, 1, l = 20),
#'      probability = TRUE, main = "Simulated data from g_vMF", xlab = "t")
#' t <- seq(-1, 1, by = 0.01)
#' lines(t, g_tilde(t = t, p = p, kappa = kappa))
#'
#' # Cosine density as a function of the dimension
#' M <- 100
#' col <- viridisLite::viridis(M)
#' plot(t, g_tilde(t = t, p = 2, kappa = kappa), col = col[2], type = "l",
#'      ylab = "Density")
#' for (p in 3:M) {
#'   lines(t, g_tilde(t = t, p = p, kappa = kappa), col = col[p])
#' }
#' @seealso \code{\link{tangent-vMF}}.
#' @name vMF


#' @rdname vMF
#' @export
d_vMF <- function(x, mu, kappa, log = FALSE) {

  # x as matrix
  if (is.null(dim(x))) {

    x <- rbind(x)

  }

  # Detect edge cases
  x <- check_unit_norm(x = x, warnings = TRUE)
  mu <- check_unit_norm(x = mu, warnings = TRUE)

  # Dimension
  p <- ncol(x)
  if (p != length(mu)) {

    stop("x and mu do not have the same dimension.")

  }

  # Density
  log_dens <- c_vMF(p = p, kappa = kappa, log = TRUE) + kappa * x %*% mu
  return(switch(log + 1, exp(log_dens), log_dens))

}


#' @rdname vMF
#' @export
c_vMF <- function(p, kappa, log = FALSE) {

  # Check kappa
  if (any(kappa < 0)) {

    stop("kappa must be non-negative.")

  }

  # Log-constant
  log_c_vMF <- (0.5 * (p - 2)) * log(kappa) - (0.5 * p) * log(2 * pi) -
    kappa - log(besselI(nu = 0.5 * (p - 2), x = kappa, expon.scaled = TRUE))
  log_c_vMF[kappa == 0] <- -w_p(p = p, log = TRUE)
  return(switch(log + 1, exp(log_c_vMF), log_c_vMF))

}


#' @rdname vMF
#' @export
r_vMF <- function(n, mu, kappa) {

  # Detect edge cases
  mu <- check_unit_norm(x = mu, warnings = TRUE)
  if (kappa < 0) {

    stop("kappa must be non-negative.")

  }

  # Dimension
  p <- length(mu)

  # Separate the uniform and the non-uniform cases
  if (kappa == 0) {

    samp <- r_unif_sphere(n = n, p = p)

  } else {

    # p = 1 is the degenerate case
    if (p > 1) {

      r_V <- function(n) r_g_vMF(n = n, p = p, kappa = kappa)
      r_U <- function(n) r_unif_sphere(n = n, p = p - 1)
      samp <- r_tang_norm(n = n, theta = mu, r_V = r_V, r_U = r_U)

    } else {

      samp <- sample(x = c(-1, 1), size = n, replace = TRUE,
                     prob = d_vMF(x = cbind(c(-1, 1)), mu = mu, kappa = kappa))

    }

  }

  # Sample
  return(samp)

}


#' @rdname vMF
#' @export
g_vMF <- function(t, p, kappa, scaled = TRUE, log = FALSE) {

  # Check kappa
  if (kappa < 0) {

    stop("kappa must be non-negative.")

  }

  # Scaled angular function
  g_c <- ifelse(scaled, c_vMF(p = p, kappa = kappa, log = TRUE), 0) +
    kappa * t
  g_c[abs(t) > 1] <- -Inf
  return(switch(log + 1, exp(g_c), g_c))

}


#' @rdname vMF
#' @export
r_g_vMF <- function(n, p, kappa) {

  # Stop if degenerate cases that will cause overflows
  if (n < 1) {

    stop("n has to be an integer larger or equal to 1")

  }
  if (p < 1) {

    stop("p has to be an integer larger or equal to 1")

  }
  if (kappa < 0) {

    stop("kappa has to be non-negative")

  }

  # Call C++ function
  return(.r_g_vMF_Cpp(n = n, p = p, kappa = kappa))

}
