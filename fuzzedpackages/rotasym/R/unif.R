

#' @title Uniform distribution on the hypersphere
#'
#' @description Density and simulation of the uniform distribution on
#' \eqn{S^{p-1}:=\{\mathbf{x}\in R^p:||\mathbf{x}||=1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p\ge 1}. The density is just the
#' inverse of the surface area of \eqn{S^{p-1}}, given by
#' \deqn{\omega_p:=2\pi^{p/2}/\Gamma(p/2).}{
#' \omega_p := 2\pi^{p/2} / \Gamma(p/2).}
#'
#' @param x locations in \eqn{S^{p-1}} to evaluate the density. Either a
#' matrix of size \code{c(nx, p)} or a vector of length \code{p}. Normalized
#' internally if required (with a \code{warning} message).
#' @param n sample size, a positive integer.
#' @param p dimension of the ambient space \eqn{R^p} that contains
#' \eqn{S^{p-1}}. A positive integer.
#' @param log flag to indicate if the logarithm of the density (or the
#' normalizing constant) is to be computed.
#' @return
#' Depending on the function:
#' \itemize{
#' \item \code{d_unif_sphere}: a vector of length \code{nx} or \code{1} with
#' the evaluated density at \code{x}.
#' \item \code{r_unif_sphere}: a matrix of size \code{c(n, p)} with the
#' random sample.
#' \item \code{w_p}: the surface area of \eqn{S^{p-1}}.
#' }
#' @details
#' If \eqn{p = 1}, then \eqn{S^{0} = \{-1, 1\}} and the "surface area" is
#' \eqn{2}. The function \code{w_p} is vectorized on \code{p}.
#' @author Eduardo García-Portugués, Davy Paindaveine, and Thomas Verdebout.
#' @examples
#' ## Area of S^{p - 1}
#'
#' # Areas of S^0, S^1, and S^2
#' w_p(p = 1:3)
#'
#' # Area as a function of p
#' p <- 1:20
#' plot(p, w_p(p = p), type = "o", pch = 16, xlab = "p", ylab = "Area",
#'      main = expression("Surface area of " * S^{p - 1}), axes = FALSE)
#' box()
#' axis(1, at = p)
#' axis(2, at = seq(0, 34, by = 2))
#'
#' ## Simulation and density evaluation for p = 1, 2, 3
#'
#' # p = 1
#' n <- 500
#' x <- r_unif_sphere(n = n, p = 1)
#' barplot(table(x) / n)
#' head(d_unif_sphere(x))
#'
#' # p = 2
#' x <- r_unif_sphere(n = n, p = 3)
#' plot(x)
#' head(d_unif_sphere(x))
#'
#' # p = 3
#' x <- r_unif_sphere(n = n, p = 3)
#' rgl::plot3d(x)
#' head(d_unif_sphere(x))
#' @name unif


#' @rdname unif
#' @export
d_unif_sphere <- function(x, log = FALSE) {

  # x as matrix
  if (is.null(dim(x))) {

    x <- rbind(x)

  }

  # Dimension
  p <- ncol(x)

  # Detect edge case x = 0 and renormalize if necessary
  x <- check_unit_norm(x = x, warnings = TRUE)

  # Density
  dens <- rep(NA, nrow(x))
  log_dens <- -w_p(p = p, log = TRUE)
  dens[complete.cases(x)] <- ifelse(log, log_dens, exp(log_dens))
  return(dens)

}


#' @rdname unif
#' @export
r_unif_sphere <- function(n, p) {

  if (p == 1) {

    # Sample {-1, 1} uniformly
    return(cbind(sample(x = c(-1, 1), size = n, replace = TRUE)))

  } else {

    # Project a N_p(0, I_p) to the sphere
    x <- matrix(rnorm(n = n * p), nrow = n, ncol = p, byrow = TRUE)
    return(x / sqrt(rowSums(x * x)))

  }

}


#' @rdname unif
#' @export
w_p <- function(p, log = FALSE) {

  # Check dimension
  if (any(p < 1)) {

    stop("p must be an integer larger or equal to 1.")

  }

  # Log-surface area
  log_w_p <- log(2) + 0.5 * p * log(pi) - lgamma(0.5 * p)

  # Exponentiated (avoid ifelse as it does not play well with vector results)
  return(switch(log + 1, exp(log_w_p), log_w_p))

}
