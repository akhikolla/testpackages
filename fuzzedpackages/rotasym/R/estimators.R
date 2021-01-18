

#' @title Estimators for the axis of rotational symmetry
#' \eqn{\boldsymbol\theta}{\theta}
#'
#' @description Estimation of the axis of rotational symmetry
#' \eqn{\boldsymbol{\theta}}{\theta} of a rotational symmetric unit-norm
#' random vector \eqn{\mathbf{X}}{X} in
#' \eqn{S^{p-1}:=\{\mathbf{x}\in R^p:||\mathbf{x}||=1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p \ge 2}, from a
#' hyperspherical sample \eqn{\mathbf{X}_1,\ldots,\mathbf{X}_n\in S^{p-1}}{
#' X_1, \ldots, X_n \in S^{p-1}}.
#' @inheritParams test_rotasym
#' @return A vector of length \code{p} with an estimate for
#' \eqn{\boldsymbol{\theta}}{\theta}.
#' @details
#' The \code{spherical_mean} estimator computes the sample mean of
#' \eqn{\mathbf{X}_1,\ldots,\mathbf{X}_n}{X_1, \ldots, X_n} and normalizes it
#' by its norm (if the norm is different from zero). It estimates consistently
#' \eqn{\boldsymbol{\theta}}{\theta} for rotational symmetric models based on
#' \link[=tang-norm-decomp]{angular functions} \eqn{g} that are
#' monotone increasing.
#'
#' The estimator in \code{spherical_loc_PCA} is based on the fact that, under
#' rotational symmetry, the expectation of
#' \eqn{\mathbf{X}\mathbf{X}'}{XX'} is
#' \eqn{a\boldsymbol{\theta}\boldsymbol{\theta}' +
#' b(\mathbf{I}_p - \boldsymbol{\theta}\boldsymbol{\theta}')}{
#' a \theta\theta'+ b (I_p - \theta\theta')}
#' for certain constants \eqn{a,b \ge 0}. Therefore,
#' \eqn{\boldsymbol{\theta}}{\theta} is the eigenvector with unique
#' multiplicity of the expectation of \eqn{\mathbf{X}\mathbf{X}'}{XX'}. Its
#' use is recommended if the rotationally symmetric data is not unimodal.
#' @author Eduardo García-Portugués, Davy Paindaveine, and Thomas Verdebout.
#' @references
#' García-Portugués, E., Paindaveine, D., Verdebout, T. (2020) On optimal tests
#' for rotational symmetry against new classes of hyperspherical distributions.
#' \emph{Journal of the American Statistical Association}, to appear.
#' \url{https://doi.org/10.1080/01621459.2019.1665527}
#' @examples
#' # Sample from a vMF
#' n <- 200
#' p <- 10
#' theta <- c(1, rep(0, p - 1))
#' set.seed(123456789)
#' data <- r_vMF(n = n, mu = theta, kappa = 3)
#' theta_mean <- spherical_mean(data)
#' theta_PCA <- spherical_loc_PCA(data)
#' sqrt(sum((theta - theta_mean)^2)) # More efficient
#' sqrt(sum((theta - theta_PCA)^2))
#'
#' # Sample from a mixture of antipodal vMF's
#' n <- 200
#' p <- 10
#' theta <- c(1, rep(0, p - 1))
#' set.seed(123456789)
#' data <- rbind(r_vMF(n = n, mu = theta, kappa = 3),
#'               r_vMF(n = n, mu = -theta, kappa = 3))
#' theta_mean <- spherical_mean(data)
#' theta_PCA <- spherical_loc_PCA(data)
#' sqrt(sum((theta - theta_mean)^2))
#' sqrt(sum((theta - theta_PCA)^2)) # Better suited in this case
#' @name estimators


#' @rdname estimators
#' @export
spherical_mean <- function(data) {

  # Check if data is spherical
  data <- check_unit_norm(x = data, warnings = TRUE)

  # Standardize spherical mean
  mu <- colMeans(data, na.rm = TRUE)
  norm <- sqrt(sum(mu * mu))
  if (norm < sqrt(.Machine$double.eps)) {

    warning("The directional mean is exactly zero")

  } else {

    mu <- mu / norm

  }
  return(mu)

}


#' @rdname estimators
#' @export
spherical_loc_PCA <- function(data) {

  # Remove NA's
  data <- data[complete.cases(data), ]

  # Check if data is spherical
  data <- check_unit_norm(x = data, warnings = TRUE)

  # Sample size and dimension
  n <- nrow(data)
  p <- ncol(data)

  # Eigen decomposition of the covariance matrix
  eig <- eigen(crossprod(data) / n, symmetric = TRUE)

  # Which one is the "outlier" (has multiplicity one)
  ind <- which.max(abs(c(
    eig$values[1] - mean(eig$values[-1]),
    eig$values[p] - mean(eig$values[-p])
  )))

  # Select the corresponding eigenvector
  theta <- switch(ind, "1" = eig$vectors[, 1], "2" = eig$vectors[, p])
  return(theta)

}
