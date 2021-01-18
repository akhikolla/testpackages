#' Informative GEV prior on a probability scale
#'
#' Constructs an informative prior for GEV parameters (\eqn{\mu, \sigma, \xi}),
#' constructed on the probability scale. For information about how to set this
#' prior see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param quant A numeric vector of length 3 containing quantiles
#'   (\eqn{q_1, q_2, q_3}{q1, q2, q3}) such that
#'   \eqn{q_1 < q_2 < q_3}{q1 < q2 < q3}. If the values
#'   in \code{quant} are not ordered from smallest to largest then they
#'   will be ordered inside \code{set_prior} without warning.
#' @param alpha A numeric vector of length 4.  Parameters specifying a
#'   prior distribution for probabilities related to the quantiles in
#'   \code{quant}.  See \strong{Details} below.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @details A prior for GEV parameters \eqn{(\mu, \sigma, \xi)},
#'   based on Crowder (1992).  This construction is typically used to set
#'   an informative prior, based on specified quantiles
#'   \eqn{q_1, q_2, q_3}{q1, q2, q3}.
#'   There are two interpretations of the parameter vector
#'   \code{alpha} = \eqn{(\alpha_1, \alpha_2, \alpha_3, \alpha_4)}:
#'   as the parameters of beta distributions for ratio of exceedance
#'   probabilities
#'   \href{https://doi.org/10.1201/b19721}{(Stephenson, 2016)}
#'   and as the parameters of
#'   a Dirichlet distribution for differences between non-exceedance
#'   probabilities
#'   \href{https://doi.org/10.1111/rssc.12159}{(Northrop et al., 2017)}.
#'   See these publications for details.
#' @return The log of the prior density.
#' @seealso \code{\link{set_prior}} for setting a prior distribution.
#' @seealso \code{\link{rpost}} and \code{\link{rpost_rcpp}} for sampling
#'   from an extreme value posterior distribution.
#' @seealso Sets the same prior as the function
#'   \code{\link[evdbayes:prior]{prior.prob}} in the evdbayes package.
#' @references Crowder, M. (1992) Bayesian priors based on parameter
#'   transformation using the distribution function
#'   \emph{Ann. Inst. Statist. Math.}, \strong{44}, 405-416.
#'   \url{https://link.springer.com/article/10.1007/BF00050695}.
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{https://doi.org/10.1111/rssc.12159}
#' @references Stephenson, A. (2016) Bayesian inference for extreme value
#'   modelling.  In \emph{Extreme Value Modeling and Risk Analysis: Methods
#'   and Applications} (eds D. K. Dey and J. Yan), 257-280, Chapman and Hall,
#'   London. \url{https://doi.org/10.1201/b19721}.
#' @export
gev_prob <- function(pars, quant, alpha, min_xi = -Inf, max_xi = Inf,
                     trendsd = 0) {
  # Extract GEV parameter values.
  mu <- pars[1]
  sigma <- pars[2]
  xi <- pars[3]
  # prior is zero if scale parameter is non-positive
  if (sigma <= 0) {
    return(-Inf)
  }
  y <- (quant - mu) / sigma
  lin <- 1 + xi * y
  # prior is zero if any component of lin is not positive.
  if (any(lin <= 0)) {
    return(-Inf)
  }
  #
  # If abs(xi) < xi_tol then xi is treated as being close to zero.
  xi_tol <- 1.e-6
  if (abs(xi) > xi_tol) {
    h <- y / xi - lin * log(lin) / xi ^ 2
  } else {
    j <- 0:4
    h_fun <- function(x) {
      sum((-1) ^ (j+1) * x ^ (j + 2) * xi ^ j / (j + 1) / (j + 2))
    }
    h <- sapply(y, h_fun)
  }
  #
  mat <- matrix(1, 3, 3)
  mat[2, ] <- y
  mat[3, ] <- h
  log_det_mat <- determinant(mat, logarithm = TRUE)$modulus
  # log-density of GEV at q
  log_g <- dgev(x = quant, loc = mu, scale = sigma, shape = xi, log = TRUE)
  pq <- pgev(q = quant, loc = mu, scale = sigma, shape = xi)
  # log-Jacobian
  log_J <- log(sigma) + sum(log_g) + log_det_mat
  # Calculate the Dirichlet log-prior
  # add the endpoints of the probability scale
  pq <- c(0, pq, 1)
  # differences between the pq values
  diff_pq <- diff(pq)
  if (any(diff_pq <= 0)) {
    return(-Inf)
  }
  log_dir_prior <- sum((alpha - 1) * log(diff_pq))
  # Combine the Jacobian with the Dirichlet prior
  val <- log_J + log_dir_prior
  return(val)
}

#' Informative GEV prior on a quantile scale
#'
#' Informative GEV prior for GEV parameters (\eqn{\mu, \sigma, \xi})
#' constructed on the quantile scale.  For information about how to set this
#' prior see \code{\link{set_prior}}.
#'
#' @param pars A numeric vector of length 3.
#'   GEV parameters (\eqn{\mu, \sigma, \xi}).
#' @param prob A numeric vector of length 3 containing exceedance
#'   probabilities (\eqn{p_1, p_2, p_3}{p1, p2, p3}) such that
#'   \eqn{p_1 > p_2 > p_3}{p1 > p2 > p3}.
#'   If the values in \code{quant} are not ordered from largest to smallest
#'   then they will be ordered inside \code{set_prior} without warning.
#' @param shape,scale Numeric vectors of length 3. Shape and scale
#'   parameters specifying (independent) gamma prior distributions placed
#'   on the differences between the quantiles corresponding to the
#'   probabilities given in \code{prob}.
#' @param min_xi  A numeric scalar.  Prior lower bound on \eqn{\xi}.
#' @param max_xi  A numeric scalar.  Prior upper bound on \eqn{\xi}.
#' @param trendsd  Has no function other than to achieve compatibility with
#'   function in the evdbayes package.
#' @return The log of the prior density.
#' @details See Coles and Tawn (1996) and/or Stephenson (2016) for details.
#'
#'   Note that the lower end point of the distribution of the distribution
#'   of the variable in question is assumed to be equal to zero.
#'   If this is not the case then the user should shift the data to
#'   ensure that this is true.
#' @references Coles, S. G. and Tawn, J. A. (1996) A Bayesian analysis of
#'   extreme rainfall data. \emph{Appl. Statist.}, \strong{45}, 463-478.
#' @references Stephenson, A. (2016) Bayesian inference for extreme value
#'   modelling.  In \emph{Extreme Value Modeling and Risk Analysis: Methods
#'   and Applications} (eds D. K. Dey and J. Yan), 257-280, Chapman and Hall,
#'   London. \url{https://doi.org/10.1201/b19721}.
#' @export
gev_quant <- function(pars, prob, shape, scale, min_xi = -Inf, max_xi = Inf,
                      trendsd = 0){
  # Extract GEV parameter values.
  mu <- pars[1]
  sigma <- pars[2]
  xi <- pars[3]
  # prior is zero if scale parameter is non-positive
  if (sigma <= 0) {
    return(-Inf)
  }
  # Calculate quantiles. Note: prob contains exceedance probabilities
  quant <- qgev(p = 1 - prob, loc = mu, scale = sigma, shape = xi)
  # If the combination of (mu, sigma, xi) is such that any of the quantiles
  # (quant[1] is the smallest) are non-positive then return -Inf.
  if (quant[1] <= 0) {
    return(-Inf)
  }
  y <- (quant - mu) / sigma
  lin <- 1 + xi * y
  # prior is zero if any component of lin is not positive.
  if (any(lin <= 0)) {
    return(-Inf)
  }
  #
  # If abs(xi) < xi_tol then xi is treated as being close to zero.
  xi_tol <- 1.e-6
  if (abs(xi) > xi_tol) {
    h <- y / xi - lin * log(lin) / xi ^ 2
  } else {
    j <- 0:4
    h_fun <- function(x) {
      sum((-1) ^ (j+1) * x ^ (j + 2) * xi ^ j / (j + 1) / (j + 2))
    }
    h <- sapply(y, h_fun)
  }
  #
  mat <- matrix(1, 3, 3)
  mat[2, ] <- y
  mat[3, ] <- h
  log_det_mat <- determinant(mat, logarithm = TRUE)$modulus
  # log-Jacobian
  log_J <- log(sigma) + log_det_mat
  # Calculate the gamma log-prior
  # differences between the quant values
  diff_quant <- diff(c(0, quant))
  log_gamma_prior <- sum((shape - 1) * log(diff_quant) - diff_quant / scale)
  # Combine the Jacobian with the Dirichlet prior
  val <- log_J + log_gamma_prior
  return(val)
}
