

#' @title Angular central Gaussian distribution
#'
#' @description Density and simulation of the Angular Central Gaussian (ACG)
#' distribution on
#' \eqn{S^{p-1}:=\{\mathbf{x}\in R^p:||\mathbf{x}||=1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p\ge 1}. The density at
#' \eqn{\mathbf{x} \in S^{p-1}}{x \in S^{p-1}}, \eqn{p\ge 2}, is given by
#' \deqn{c^{\mathrm{ACG}}_{p,\boldsymbol{\Lambda}}
#' (\mathbf{x}' \boldsymbol{\Lambda}^{-1} \mathbf{x})^{-p/2}
#' \quad\mathrm{with}\quad c^{\mathrm{ACG}}_{p,\boldsymbol{\Lambda}}:=
#' 1 / (\omega_p |\boldsymbol{\Lambda}|^{1/2})}{
#' c^{ACG}_{p,\Lambda} (x' \Lambda^{-1} x)^{-p/2} with
#' c^{ACG}_{p,\Lambda} := 1 / (\omega_p |\Lambda|^{1/2}).}
#' where \eqn{\boldsymbol{\Lambda}}{\Lambda} is the shape matrix, a
#' \eqn{p\times p}{p x p} symmetric and positive definite matrix, and
#' \eqn{\omega_p} is the surface area of \eqn{S^{p-1}}.
#'
#' @inheritParams unif
#' @param Lambda the shape matrix \eqn{\boldsymbol{\Lambda}}{\Lambda} of the
#' ACG. A symmetric and positive definite matrix of size \code{c(p, p)}.
#' @return
#' Depending on the function:
#' \itemize{
#' \item \code{d_ACG}: a vector of length \code{nx} or \code{1} with the
#' evaluated density at \code{x}.
#' \item \code{r_ACG}: a matrix of size \code{c(n, p)} with the random sample.
#' \item \code{c_ACG}: the normalizing constant.
#' }
#' @details
#' Due to the projection of the ACG, the shape matrix
#' \eqn{\boldsymbol{\Lambda}}{\Lambda} is only identified up to a constant,
#' that is, \eqn{\boldsymbol{\Lambda}}{\Lambda} and
#' \eqn{c\boldsymbol{\Lambda}}{c\Lambda} give the same ACG distribution.
#' Usually, \eqn{\boldsymbol{\Lambda}}{\Lambda} is normalized to have trace
#' equal to \eqn{p}.
#'
#' \code{c_ACG} is vectorized on \code{p}. If \eqn{p = 1}, then the ACG is the
#' uniform distribution in the set \eqn{\{-1, 1\}}{{-1, 1}}.
#' @author Eduardo García-Portugués, Davy Paindaveine, and Thomas Verdebout.
#' @references
#' Tyler, D. E. (1987). Statistical analysis for the angular central Gaussian
#' distribution on the sphere. \emph{Biometrika}, 74(3):579--589.
#' \url{https://doi.org/10.1093/biomet/74.3.579}
#' @examples
#' # Simulation and density evaluation for p = 2
#' Lambda <- diag(c(5, 1))
#' n <- 1e3
#' x <- r_ACG(n = n, Lambda = Lambda)
#' col <- viridisLite::viridis(n)
#' r <- runif(n, 0.95, 1.05) # Radius perturbation to improve visualization
#' plot(r * x, pch = 16, col = col[rank(d_ACG(x = x, Lambda = Lambda))])
#'
#' # Simulation and density evaluation for p = 3
#' Lambda <- rbind(c(5, 1, 0.5),
#'                 c(1, 2, 1),
#'                 c(0.5, 1, 1))
#' x <- r_ACG(n = n, Lambda = Lambda)
#' rgl::plot3d(x, col = col[rank(d_ACG(x = x, Lambda = Lambda))], size = 5)
#' @seealso \code{\link{tangent-elliptical}}.
#' @name ACG


#' @rdname ACG
#' @export
d_ACG <- function(x, Lambda, log = FALSE) {

  # x as matrix
  if (is.null(dim(x))) {

    x <- rbind(x)

  }

  # Dimension
  p <- ncol(x)
  if (p != sqrt(length(Lambda))) {

    stop("x and Lambda do not have the same dimension.")

  }
  if (p == 1) {

    # Log-density
    log_dens <- d_unif_sphere(x = x, log = TRUE)

  } else {

    # Detect edge case x = 0 and renormalize if necessary
    x <- check_unit_norm(x = x, warnings = TRUE)

    # Log-density
    log_dens <- c_ACG(p = p, Lambda = Lambda, log = TRUE) -
      0.5 * (p - 1) * log(rowSums((x %*% solve(Lambda)) * x))

  }
  return(switch(log + 1, exp(log_dens), log_dens))

}


#' @rdname ACG
#' @export
c_ACG <- function(p, Lambda, log = FALSE) {

  # Check symmetry of Lambda
  if (!isSymmetric(Lambda, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {

    stop("Lambda must be a symmetric matrix")

  }

  # Log-determinant using Cholesky decomposition (generates an error if Lambda
  # is not positive definite)
  log_det <- 2 * sum(log(diag(chol(Lambda))))
  log_c_ACG <- - (w_p(p = p, log = TRUE) + 0.5 * log_det)
  return(switch(log + 1, exp(log_c_ACG), log_c_ACG))

}


#' @rdname ACG
#' @export
r_ACG <- function(n, Lambda) {

  # Dimension
  p <- sqrt(length(Lambda))

  # Simulation (chol() generates an error if Lambda is not positive definite)
  x <- matrix(rnorm(n = n * p), nrow = n, ncol = p, byrow = TRUE) %*%
    chol(Lambda)

  # Projection
  return(x / sqrt(rowSums(x * x)))

}
