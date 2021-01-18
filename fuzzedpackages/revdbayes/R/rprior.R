# ================================ rprior_quant ===============================

#' Prior simulation of GEV parameters - prior on quantile scale
#'
#' Simulates from the prior distribution for GEV parameters proposed in
#' Coles and Tawn (1996), based on independent gamma priors for differences
#' between quantiles.
#'
#' @param n A numeric scalar. The size of sample required.
#' @param prob A numeric vector of length 3. Exceedance probabilities
#'   corresponding to the quantiles used to specify the prior distribution.
#'   The values should \emph{decrease} with the index of the vector.
#'   If not, the values in \code{prob} will be sorted into decreasing order
#'   without warning.
#' @param shape A numeric vector of length 3. Respective shape parameters of
#'   the gamma priors for the quantile differences.
#' @param scale A numeric vector of length 3. Respective scale parameters of
#'   the gamma priors for the quantile differences.
#' @param lb A numeric scalar.  If this is not \code{NULL} then the simulation
#'   is constrained so that \code{lb} is an approximate lower bound on the
#'   GEV variable.  Specifically, only simulated GEV parameter values for
#'   which the 100\code{lb_prob}\% quantile is greater than \code{lb} are
#'   retained.
#' @param lb_prob A numeric scalar.  The non-exceedance probability involved
#'   in the specification of \code{lb}.  Must be in (0,1).  If \code{lb=NULL}
#'   then \code{lb_prob} is not used.
#' @details The simulation is based on the way that the prior is constructed.
#'   See Coles and Tawn (1996), Stephenson (2016) or the evdbayes user guide
#'   for details of the construction of the prior. First, the quantile
#'   differences are simulated from the specified gamma distributions.
#'   Then the simulated quantiles are calculated. Then the GEV location,
#'   scale and shape parameters that give these quantile values are found,
#'   by solving numerically a set of three non-linear equations in which the
#'   GEV quantile function evaluated at the values in \code{prob} is equated
#'   to the simulated quantiles.  This is reduced to a one-dimensional
#'   optimisation over the GEV shape parameter.
#' @return An \code{n} by 3 numeric matrix.
#' @seealso \code{\link[evdbayes:prior]{prior.quant}} to set this prior
#'   using the evdbayes package.
#' @seealso \code{\link[evdbayes:prior]{posterior}}: evdbayes function that
#'   can sample from this prior distribution (using MCMC) if the argument
#'   \code{lh = "none"} is given.
#' @seealso \code{\link{rpost}} and \code{\link{rpost_rcpp}} for sampling
#'   from an extreme value posterior distribution.
#' @references Coles, S. G. and Tawn, J. A. (1996) A Bayesian analysis of
#'   extreme rainfall data. \emph{Appl. Statist.}, \strong{45}, 463-478.
#' @references Stephenson, A. 2016. Bayesian Inference for Extreme Value
#'   Modelling. In \emph{Extreme Value Modeling and Risk Analysis: Methods and
#'   Applications}, edited by D. K. Dey and J. Yan, 257-80. London:
#'   Chapman and Hall. \url{https://doi.org/10.1201/b19721}
#' @examples
#' pr <- 10 ^ -(1:3)
#' sh <- c(38.9, 7.1, 47)
#' sc <- c(1.5, 6.3, 2.6)
#' x <- rprior_quant(n = 1000, prob = pr, shape = sh, scale = sc)
#' x <- rprior_quant(n = 1000, prob = pr, shape = sh, scale = sc, lb = 0)
#' @export
rprior_quant <- function(n, prob, shape, scale, lb = NULL, lb_prob = 0.001){
  #
  pq <- sort(prob, decreasing = TRUE)
  # n_sim is the number of (extra) values that need to be simulated.
  n_sim <- n
  # Matrix in which to store the results.
  pars <- matrix(NA, nrow = n, ncol = 3)
  #
  next_row <- 1
  while (n_sim > 0) {
    #
    # Simulate qp_tildes (differences between quantiles)
    #
    qp1_tilde <- stats::rgamma(n, shape = shape[1], scale = scale[1])
    qp2_tilde <- stats::rgamma(n, shape = shape[2], scale = scale[2])
    qp3_tilde <- stats::rgamma(n, shape = shape[3], scale = scale[3])
    #
    # Transform to qps
    #
    qp1 <- qp1_tilde
    qp2 <- qp1_tilde + qp2_tilde
    qp3 <- qp1_tilde + qp2_tilde + qp3_tilde
    qp_sim <- cbind(qp1, qp2, qp3)
    #
    # Transform to (mu, sigma, xi)
    #
    temp <- t(apply(qp_sim, 1, quantile_to_gev, prob = prob))
    #
    # If lb is not NULL then check for admissibility of the simulated values.
    #
    if (!is.null(lb)) {
      check <- qgev(p = lb_prob, loc = temp[, 1], scale = temp[, 2],
                    shape = temp[, 3])
      # Retain only the values with 100lb_prob% quantile > lb.
      temp <- temp[check > lb, , drop = FALSE]
      # ... up to a maximum of the remaining number that we require.
      n_ret <- min(nrow(temp), n_sim)
      # Add values to pars.
      pars[next_row:(next_row + n_ret - 1), ] <- temp[1:n_ret, ]
      next_row <- next_row + n_ret
      # Reduce the number that we require.
      n_sim <- n_sim - n_ret
    } else {
      pars <- temp
      n_sim <- 0
    }
  }
  return(pars)
}

# ================================ rprior_prob ===============================

#' Prior simulation of GEV parameters - prior on probability scale
#'
#' Simulates from the prior distribution for GEV parameters based on
#' Crowder (1992), in which independent beta priors are specified
#' for ratios of probabilities (which is equivalent to
#' a Dirichlet prior on differences between these probabilities).
#'
#' @param n A numeric scalar. The size of sample required.
#' @param quant A numeric vector of length 3.  Contains quantiles
#'   \eqn{q_1, q_2, q_3}.  A prior distribution is placed on the
#'   non-exceedance (\code{exc = FALSE}) or exceedance (\code{exc = TRUE})
#'   probabilities corresponding to these quantiles.
#'   The values should \emph{increase} with the index of the vector.
#'   If not, the values in \code{quant} will be sorted into increasing order
#'   without warning.
#' @param alpha A numeric vector of length 4. Parameters of the Dirichlet
#'   distribution for the exceedance probabilities.
#' @param exc A logical scalar.  Let \eqn{M} be the GEV variable,
#'   \eqn{r_q = P(M \leq q)}{r_q = P(M <= q)},
#'   \eqn{p_q = P(M > q) = 1 - r_q} and
#'   \code{quant} = (\eqn{q_1, q_2, q_3}).
#'   If \code{exc = FALSE} then a Dirichlet(\code{alpha}) distribution is
#'   placed on
#'   \eqn{(r_{q_1}, r_{q_2} - r_{q_1}, r_{q_3} - r_{q_2}, 1 - r_{q_3})}{%
#'         (r_q1, r_q2 - r_q1, r_q3 - r_q2, 1 - r_q3)}, as in
#'   Northrop et al. (2017).
#'   If \code{exc = TRUE} then a Dirichlet(\code{alpha}) distribution
#'   is placed on
#'   \eqn{(1 - p_{q_1}, p_{q_1} - p_{q_2}, p_{q_2} - p_{q_3}, p_{q_3})}{%
#'        (1 - p_q1, p_q1 - p_q2, p_q2 - p_q3, p_q3)}, where
#'   \eqn{p_q = P(M > q)}, as in Stephenson (2016).
#' @param lb A numeric scalar.  If this is not \code{NULL} then the simulation
#'   is constrained so that \code{lb} is an approximate lower bound on the
#'   GEV variable.  Specifically, only simulated GEV parameter values for
#'   which the 100\code{lb_prob}\% quantile is greater than \code{lb} are
#'   retained.
#' @param lb_prob A numeric scalar.  The non-exceedance probability involved
#'   in the specification of \code{lb}.  Must be in (0,1).  If \code{lb=NULL}
#'   then \code{lb_prob} is not used.
#' @details The simulation is based on the way that the prior is constructed.
#'   See
#'   \href{https://doi.org/10.1201/b19721}{Stephenson (1996)}
#'   the evdbayes user guide or Northrop et al. (2017)
#'   \href{https://doi.org/10.1111/rssc.12159}{Northrop et al. (2017)}
#'   for details of the construction of the prior.  First, differences between
#'   probabilities are simulated from a Dirichlet distribution. Then the GEV
#'   location, scale and shape parameters that correspond to these quantile
#'   values are found, by solving numerically a set of three non-linear
#'   equations in which the GEV quantile function evaluated at the simulated
#'   probabilities is equated to the quantiles in \code{quant}.
#'   This is reduced to a one-dimensional optimisation over the GEV shape
#'   parameter.
#' @return An \code{n} by 3 numeric matrix.
#' @seealso \code{\link[evdbayes:prior]{prior.prob}} to set this prior using
#'   the evdbayes package.
#' @seealso \code{\link[evdbayes]{posterior}}: evdbayes function that can
#'   sample from this prior distribution (using MCMC) if the argument
#'   \code{lh = "none"} is given.
#' @seealso \code{\link{rpost}} and \code{\link{rpost_rcpp}} for sampling
#'   from an extreme value posterior distribution.
#' @references Crowder, M. (1992) Bayesian priors based on parameter
#'   transformation using the distribution function.
#'   \emph{Ann. Inst. Statist. Math.}, \strong{44}(3), 405-416.
#'   \url{https://link.springer.com/article/10.1007/BF00050695}
#' @references Stephenson, A. 2016. Bayesian Inference for Extreme Value
#'   Modelling. In \emph{Extreme Value Modeling and Risk Analysis: Methods and
#'   Applications}, edited by D. K. Dey and J. Yan, 257-80. London:
#'   Chapman and Hall. \url{https://doi.org/10.1201/b19721}
#' @references Northrop, P. J., Attalides, N. and Jonathan, P. (2017)
#'   Cross-validatory extreme value threshold selection and uncertainty
#'   with application to ocean storm severity.
#'   \emph{Journal of the Royal Statistical Society Series C: Applied
#'   Statistics}, \strong{66}(1), 93-120.
#'   \url{https://doi.org/10.1111/rssc.12159}
#' @examples
#' quant <- c(85, 88, 95)
#' alpha <- c(4, 2.5, 2.25, 0.25)
#' x <- rprior_prob(n = 1000, quant = quant, alpha = alpha, exc = TRUE)
#' x <- rprior_prob(n = 1000, quant = quant, alpha = alpha, exc = TRUE, lb = 0)
#' @export
rprior_prob <- function(n, quant, alpha, exc = FALSE, lb = NULL,
                        lb_prob = 0.001){
  #
  quant <- sort(quant, decreasing = FALSE)
  # n_sim is the number of (extra) values that need to be simulated.
  n_sim <- n
  # Matrix in which to store the results.
  pars <- matrix(NA, nrow = n, ncol = 3)
  #
  next_row <- 1
  while (n_sim > 0) {
    #
    # Simulate differences between probabilities from a Dirichlet disribution.
    #
    if (exc) {
      pq_diff <- rDir(n = n_sim, alpha = alpha)
      # Dirichlet prior is on (1 - p_q1, p_q1 - p_q2, p_q2 - p_q3, p_q3)
      pq3 <- pq_diff[, 4]
      pq1 <- 1 - pq_diff[, 1]
      pq2 <- pq_diff[, 3] + pq3
      prob <- cbind(pq1, pq2, pq3)
    }
    else {
      rq_diff <- rDir(n = n_sim, alpha = alpha)
      # Dirichlet prior is on (r_q1, r_q2 - r_q1, r_q3 - r_q2, 1 - r_q3)
      rq1 <- rq_diff[, 1]
      rq3 <- 1 - rq_diff[, 4]
      rq2 <- rq_diff[, 2] + rq1
      # quantile_to_gev() uses exceedance probabilities.
      prob <- 1 - cbind(rq1, rq2, rq3)
    }
    #
    # Transform to (mu, sigma, xi)
    #
    temp <- t(apply(prob, 1, quantile_to_gev, quant = quant))
    #
    # If lb is not NULL then check for admissibility of the simulated values.
    #
    if (!is.null(lb)) {
      check <- qgev(p = lb_prob, loc = temp[, 1], scale = temp[, 2],
                    shape = temp[, 3])
      # Retain only the values with 100lb_prob% quantile > lb.
      temp <- temp[check > lb, , drop = FALSE]
      # ... up to a maximum of the remaining number that we require.
      n_ret <- min(nrow(temp), n_sim)
      # Add values to pars.
      pars[next_row:(next_row + n_ret - 1), ] <- temp[1:n_ret, ]
      next_row <- next_row + n_ret
      # Reduce the number that we require.
      n_sim <- n_sim - n_ret
    } else {
      pars <- temp
      n_sim <- 0
    }
  }
  return(pars)
}

# ================================ quantile_to_gev ===============================

#' Converts quantiles to GEV parameters
#'
#' Three quantiles, that is, the value of quantile and their respective
#' exceedance probabilities, are provided. This function attempts to
#' find the location, scale and shape parameters of a GEV distribution that
#' has these quantiles.
#'
#' @param quant A numeric vector of length 3. Values of the quantiles.
#'   The values should \emph{increase} with the index of the vector.
#'   If not, the values in \code{quant} will be sorted into increasing order
#'   without warning.
#' @param prob A numeric vector of length 3. Exceedance probabilities
#'   corresponding to the quantiles in \code{quant}.
#'   The values should \emph{decrease} with the index of the vector.
#'   If not, the values in \code{prob} will be sorted into decreasing order
#'   without warning.
#' @details Suppose that \eqn{G(x)} is the distribution function of
#'   a GEV(\eqn{\mu, \sigma, \xi}) distribution.  This function attempts to
#'   solve numerically the set of three non-linear equations
#'   \deqn{G(q_i) = 1 - p_i, i = 1, 2, 3}
#'   where \eqn{q_i, i=1,2,3} are the quantiles in \code{quant} and
#'   \eqn{p_i, i=1,2,3} are the exceedance probabilities in \code{prob}.
#'   This is reduced to a one-dimensional optimisation over the GEV
#'   shape parameter.
#' @return A numeric vector of length 3 containing the GEV location, scale and
#'   shape parameters.
#' @seealso \code{\link{rprior_quant}} for simulation of GEV parameters from
#'   a prior constructed on the quantile scale.
#' @examples
#' my_q <- c(15, 20, 22.5)
#' my_p <- 1-c(0.5, 0.9, 0.5^0.01)
#' x <- quantile_to_gev(quant = my_q, prob = my_p)
#' # Check
#' qgev(p = 1 - my_p, loc = x[1], scale = x[2], shape = x[3])
#' @export
quantile_to_gev <- function(quant, prob){
  pq <- sort(prob, decreasing = TRUE)
  qp <- sort(quant, decreasing = FALSE)
  #
  # Transform to (mu, sigma, xi)
  #
  f_xi <- function(xi, pq){
    xp <- -log(1 - pq)
    ifelse(abs(xi) < 1e-6, -log(xp) + xi * log(xp) ^ 2 / 2,
           (xp ^ (-xi) - 1) / xi)
  }
  #
  obfn <- function(xi, qp, pq){
    f1 <- f_xi(xi = xi, pq = pq[1])
    f2 <- f_xi(xi = xi, pq = pq[2])
    f3 <- f_xi(xi = xi, pq = pq[3])
    sigma <- (qp[3] - qp[1])/(f3 - f1)
    mu <- qp[1] - sigma * f1
    (mu + sigma * f2 - qp[2])^2
  }
  #
  xi <- stats::nlminb(0.01, obfn, qp = qp, pq = pq)$par
  f1 <- f_xi(xi = xi, pq = pq[1])
  f3 <- f_xi(xi = xi, pq = pq[3])
  sigma <- (qp[3] - qp[1])/(f3 - f1)
  mu <- qp[1] - sigma * f1
  #
  gevpars <- c(mu, sigma, xi)
  names(gevpars) <- c("location", "scale", "shape")
  return(gevpars)
}

# ================================ rDir ===============================

#' Simulation from a Dirichlet distribution
#'
#' Simulates from a Dirichlet distribution with concentration parameter
#' vector \eqn{\alpha} = (\eqn{\alpha_1}, ..., \eqn{\alpha_K}).
#'
#' @param n A numeric scalar. The size of sample required.
#' @param alpha A numeric vector.  Dirichlet concentration parameter.
#' @details The simulation is based on the property that if
#'   \eqn{Y_1, \ldots, Y_K}{Y_1, ..., Y_K} are independent, \eqn{Y_i} has a
#'   gamma(\eqn{\alpha_i}, 1) distribution and
#'   \eqn{S = Y_1 + \cdots + Y_k}{S = Y_1 + ... + Y_k}
#'   then \eqn{(Y_1, \ldots, Y_K) / S}{(Y_1, ..., Y_K) / S} has a
#'   Dirichlet(\eqn{\alpha_1}, ..., \eqn{\alpha_K}) distribution.
#'
#'   See
#' \url{https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution}
#' @return An \code{n} by \code{length(alpha)} numeric matrix.
#' @seealso \code{\link{rprior_prob}} for prior simulation of
#'   GEV parameters - prior on probability scale.
#' @references Kotz, S., Balakrishnan, N. and Johnson, N. L. (2000)
#'   \emph{Continuous Multivariate Distributions, vol. 1, Models and
#'   Applications, 2nd edn}, ch. 49. New York: Wiley.
#'   \url{https://doi.org/10.1002/0471722065}
#' @examples
#' rDir(n = 10, alpha = 1:4)
#' @export
rDir <- function(n = 1, alpha = c(1, 1)){
  y <- matrix(NA, nrow = n, ncol = length(alpha))
  for (j in 1:ncol(y)){
    y[, j] <- stats::rgamma(n, shape = alpha[j])
  }
  y / rowSums(y)
}

