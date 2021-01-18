#' Maximum likelihood estimate for \eqn{\sigma^2}, \eqn{\phi} and \eqn{\rho}.
#'
#' This function maximum likelihood estimate for \eqn{\sigma^2}, \eqn{\phi}
#'  and \eqn{\rho} in the random field model for the covariance.
#'
#' @param m_coef Matrix where each column is an observed vector
#' @param m_coord Matrix where each observation contains the latitude and
#' longitude
#' @export
#' @return Return a list with
#' \describe{
#'    \item{par}{A vector with the estimates of \eqn{\sigma^2}, \eqn{\phi} and
#'    \eqn{\rho}.}
#'    \item{m_cov}{A matrix of covariances of the estimates.}
#' }
#' @examples
#' data("datasetCanada")
#'
#' m_data <- as.matrix(datasetCanada$m_data)
#' m_coord <- as.matrix(datasetCanada$m_coord[, 1:2])
#'
#' p <- ceiling(1 + log2(nrow(m_data)))
#' m_coef <- sapply(seq_len(nrow(m_coord)), function(i) {
#'     coef_fourier(m_data[, i], p)
#' })
#' log_lik_rf(m_coef, m_coord)
log_lik_rf <- function(m_coef, m_coord) {

  m_dist <- as.matrix(dist(m_coord))

  m <- (nrow(m_coef) - 1) / 2

  m_dif <- sapply(seq_len(2 * m + 1), function(j) abs(seq_len(2 * m + 1) - j))

  fn <- function(arg) {
    as.numeric(logLikMultiNorm(m_coef, m_dist, arg[1], arg[2], arg[3]))
  }

  inicial <- c(2, 0.5, 0.25) + rnorm(3, sd = 0.05)

  fit <- nlminb(inicial, fn,
                gradient = function(x) numDeriv::grad(fn, x),
                lower = c(0.001, 0.001, -0.999), upper = c(Inf, Inf, 0.999))


  m_cov <- with(fit, solve(numDeriv::hessian(fn, par)))

  if (min(diag(m_cov)) < 0) stop("Error in computing the hessian matrix.")

  with(fit, list(par = par, m_cov = m_cov))
}


#' Kriging method for Spatial Functional Data.
#'
#' \code{geo_fkf} implements the kriging method for spatial functional datasets.
#'
#' @param m_data a tibble where each column or variable is data from a station
#' @param m_coord a tibble with two columns: latitude and longitude
#' @param new_loc a tible with one observation, where the columns or variables
#' are latitude and longitude
#' @param p order in the Fourier Polynomial
#' @param t a time series with values belonging to \eqn{[\-pi, \pi]} to
#' evaluate the estimate curve
#' @export
#'
#' @return a list with three entries: \code{estimates}, \code{Theta} and
#' \code{cov_params}
#' \describe{
#'   \item{estimates}{the estimate curve}
#'   \item{Theta}{weights (matrices) of the linear combination}
#'   \item{cov_params}{estimate \eqn{\sigma^2}, \eqn{\phi} and \eqn{\rho}}
#' }
#' @examples
#'data("datasetCanada")
#'
#'m_data <- as.matrix(datasetCanada$m_data)
#'m_coord <- as.matrix(datasetCanada$m_coord[, 1:2])
#'pos <- sample.int(nrow(m_coord), 1)
#'log_pos <- !(seq_len(nrow(m_coord)) %in% pos)
#'new_loc <- m_coord[pos, ]
#'m_coord <- m_coord[log_pos, ]
#'m_data <- m_data[, log_pos]
#'
#'geo_fkf(m_data, m_coord, new_loc)
geo_fkf <- function(m_data, m_coord, new_loc, p,
                    t = seq(from = -pi, to = pi, by = 0.01)) {

  n_loc <- nrow(m_coord)

  if (missing(p)) p <- ceiling(1 + log2(nrow(m_data)))

  m_coef <- sapply(seq_len(n_loc), function(i) {
    x <- (as.matrix(m_data))[, i]
    coef_fourier(x, p)
  })

  fit_cov <- log_lik_rf(as.matrix(m_coef), as.matrix(m_coord))
  s2 <- fit_cov$par[1]
  phi <- fit_cov$par[2]
  rho <- fit_cov$par[3]

  m_dist <- as.matrix(dist(rbind(as.matrix(m_coord),
                                 as.vector(new_loc))))
  m_dif <- sapply(seq_len(2 * p + 1),
                  function(j) abs(seq_len(2 * p + 1) - j))

  m_a <- matrix(0, nrow = (n_loc + 1) * (2 * p + 1),
                ncol = (n_loc + 1) * (2 * p + 1))

  for (i in seq_len(n_loc + 1)) {
    v_rows <- (n_loc * (2 * p + 1) + 1):((n_loc + 1) * (2 * p + 1))
    v_cols <- ((i - 1) * (2 * p + 1) + 1):(i * (2 * p + 1))
    m_a[v_rows, v_cols] <- m_a[v_cols, v_rows] <- diag(2 * p + 1)
  }

  for (i in 1:n_loc) {
    for (j in i:n_loc) {
      v_rows <- ((i - 1) * (2 * p + 1) + 1):(i * (2 * p + 1))
      v_cols <- ((j - 1) * (2 * p + 1) + 1):(j * (2 * p + 1))
      m_a[v_cols, v_rows] <-
        m_a[v_rows, v_cols] <-
        s2 * exp(-m_dist[i, j] / phi) * rho^m_dif * diag(2 * p + 1)
    }
  }

  m_b <- matrix(1, nrow = (n_loc + 1) * (2 * p + 1), ncol = 1)
  for (i in seq_len(n_loc)) {
    v_rows <- ((i - 1) * (2 * p + 1) + 1):(i * (2 * p + 1))
    m_b[v_rows] <- diag(s2 * exp(-m_dist[i, n_loc + 1] / phi) * rho^m_dif)
  }

  solution <- matrix(solve(m_a, m_b), nrow = 2 * p + 1)

  c0 <- matrix(0, nrow = nrow(m_coef), ncol = 1)
  for (i in seq_len(n_loc)) {
    c0 <- c0 + diag(solution[, i]) %*% m_coef[, i]
  }

  list(estimates = fourier_b(c0, x = t),
       Theta = solution[, seq_len(n_loc)],
       cov_params = list(s2 = s2, phi = phi, rho = rho))
}
