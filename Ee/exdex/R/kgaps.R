# =================================== kgaps ===================================
#
#' Maximum likelihood estimation for the \eqn{K}-gaps model
#'
#' Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
#' based on the \eqn{K}-gaps model for threshold inter-exceedances times of
#' Suveges and Davison (2010).
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.  See Attalides (2015) for details.
#' @details The maximum likelihood estimate of the extremal index \eqn{\theta}
#'   under the \eqn{K}-gaps model of Suveges and Davison (2010) is calculated.
#'   If \code{inc_cens = TRUE} then information from censored inter-exceedance
#'   times is included in the likelihood to be maximized, following
#'   Attalides (2015).  The form of the log-likelihood is given in the
#'   \strong{Details} section of \code{\link{kgaps_stat}}.
#'
#'   It is possible that the estimate of \eqn{\theta} is equal to 1, and also
#'   possible that it is equal to 0. \code{\link{kgaps_stat}} explains the
#'   respective properties of the data that cause these events to occur.
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{http://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @return An object (a list) of class \code{c("kgaps", "exdex")} containing
#'     \item{\code{theta} }{The maximum likelihood estimate (MLE) of
#'       \eqn{\theta}.}
#'     \item{\code{se} }{The estimated standard error of the MLE.}
#'     \item{\code{ss} }{The list of summary statistics returned from
#'       \code{\link{kgaps_stat}}.}
#'     \item{\code{k, u, inc_cens} }{The input values of \code{k},
#'       \code{u} and \code{inc_cens}.}
#'     \item{\code{call }}{The call to \code{kgaps}.}
#' @seealso \code{\link{confint.kgaps}} to estimate confidence intervals
#'   for \eqn{\theta}.
#' @seealso \code{\link{kgaps_imt}} for the information matrix test, which
#'   may be used to inform the choice of the pair (\code{u, k}).
#' @seealso \code{\link{choose_uk}} for a diagnostic plot based on
#'   \code{\link{kgaps_imt}}.
#' @seealso \code{\link{kgaps_stat}} for the calculation of sufficient
#'   statistics for the \eqn{K}-gaps model.
#' @seealso \code{\link[revdbayes]{kgaps_post}} in the
#'   \code{\link[revdbayes]{revdbayes}} package for Bayesian inference
#'   about \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{spm}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @seealso \code{\link{iwls}}: iterated weighted least squares estimator.
#' @examples
#' ### S&P 500 index
#' u <- quantile(sp500, probs = 0.60)
#' theta <- kgaps(sp500, u)
#' theta
#' summary(theta)
#'
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = 0.60)
#' theta <- kgaps(newlyn, u, k= 2)
#' theta
#' summary(theta)
#' @export
kgaps <- function(data, u, k = 1, inc_cens = FALSE) {
  Call <- match.call(expand.dots = TRUE)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (u >= max(data)) {
    stop("u must be less than max(data)")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("k must be a numeric scalar")
  }
  # Calculate sufficient statistics
  ss <- kgaps_stat(data, u, k, inc_cens)
  # If N0 = 0 then all exceedances occur singly (all K-gaps are positive)
  # and the likelihood is maximized at theta = 1.
  N0 <- ss$N0
  # If N1 = 0 then we are in the degenerate case where there is one cluster
  # (all K-gaps are zero) and the likelihood is maximized at theta = 0.
  N1 <- ss$N1
  if (N1 == 0) {
    theta_mle <- 0L
  } else if (N0 == 0) {
    theta_mle <- 1L
  } else {
    sum_qs <- ss$sum_qs
    theta_mle <- kgaps_quad_solve(N0, N1, sum_qs)
  }
  # Estimate standard error
  obs_info <- 0
  if (N0 > 0) {
    obs_info <- obs_info + N0 / (1 - theta_mle) ^ 2
  }
  if (N1 > 0) {
    obs_info <- obs_info + 2 * N1 / theta_mle ^ 2
  }
  theta_se <- sqrt(1 / obs_info)
  res <- list(theta = theta_mle, se = theta_se, ss = ss, k = k,
              u = u, inc_cens = inc_cens, call = Call)
  class(res) <- c("kgaps", "exdex")
  return(res)
}

# ================================ kgaps_stat =================================

#' Sufficient statistics for the \eqn{K}-gaps model
#'
#' Calculates sufficient statistics for the \eqn{K}-gaps model for the extremal
#' index \eqn{\theta}.
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = \max(T - K, 0)}{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.  See Attalides (2015) for details.
#' @details The sample \eqn{K}-gaps are
#'   \eqn{S_0, S_1, ..., S_{N-1}, S_N}{S_0, S_1, ..., S_(N-1), S_N},
#'   where \eqn{S_1, ..., S_{N-1}}{S_1, ..., S_(N-1)} are uncensored and
#'   \eqn{S_0} and \eqn{S_N} are censored.  Under the assumption that the
#'   \eqn{K}-gaps are independent, the log-likelihood of the \eqn{K}-gaps
#'   model is given by
#'   \deqn{l(\theta; S_0, \ldots, S_N) = N_0 \log(1 - \theta) +
#'     2 N_1 \log \theta - \theta q (S_0 + \cdots + S_N),}{%
#'     l(\theta; S_0, ..., S_N) = N_0 log(1 - \theta) + 2 N_1 log \theta -
#'     \theta q (S_0 + ... + S_N),}
#'    where \eqn{q} is the threshold exceedance probability,
#'    \eqn{N_0} is the number of sample \eqn{K}-gaps that are equal to zero and
#'    (apart from an adjustment for the contributions of \eqn{S_0} and
#'    \eqn{S_N}) \eqn{N_1} is the number of positive sample \eqn{K}-gaps.
#'    Specifically, \eqn{N_1} is equal to the number of
#'    \eqn{S_1, ..., S_{N-1}}{S_1, ..., S_(N-1)}
#'    that are positive plus \eqn{(I_0 + I_N) / 2}, where \eqn{I_0 = 1} if
#'    \eqn{S_0} is greater than zero and similarly for \eqn{I_N}.
#'    The differing treatment of uncensored and censored \eqn{K}-gaps reflects
#'    differing contributions to the likelihood.
#'    For full details see Suveges and Davison (2010) and Attalides (2015).
#'
#'    If \eqn{N_1 = 0} then we are in the degenerate case where there is one
#'    cluster (all \eqn{K}-gaps are zero) and the likelihood is maximized at
#'    \eqn{\theta = 0}.
#'
#'    If \eqn{N_0 = 0} then all exceedances occur singly (all \eqn{K}-gaps are
#'    positive) and the likelihood is maximized at \eqn{\theta = 1}.
#' @return A list containing the sufficient statistics, with components
#'     \item{\code{N0} }{the number of zero \eqn{K}-gaps}
#'     \item{\code{N1} }{contribution from non-zero \eqn{K}-gaps (see
#'       \strong{Details})}
#'     \item{\code{sum_qs} }{the sum of the (scaled) \eqn{K}-gaps, i.e.
#'       \eqn{q (S_0 + \cdots + S_N)}{q (S_0 + ... + S_N)}, where \eqn{q}
#'       is estimated by the proportion of threshold exceedances.}
#'     \item{\code{n_kgaps} }{the number of \eqn{K}-gaps, including 2
#'       censored \eqn{K}-gaps if \code{inc_cens = TRUE}.}
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#'   \url{http://discovery.ucl.ac.uk/1471121/1/Nicolas_Attalides_Thesis.pdf}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @examples
#' u <- quantile(newlyn, probs = 0.90)
#' kgaps_stat(newlyn, u)
#' @export
kgaps_stat <- function(data, u, k = 1, inc_cens = FALSE) {
  if (any(is.na(data))) {
    stop("No missing values are allowed in ''data''")
  }
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (u >= max(data)) {
    stop("u must be less than max(data)")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("k must be a numeric scalar")
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > u]
  N_u <- length(exc_u)
#  q_u <- (N_u - 1) / nx # mev
  q_u <- N_u / nx
  # Inter-exceedances times and K-gaps
  T_u <- diff(exc_u)
  S_k <- pmax(T_u - k, 0)
  # N0, N1, sum of scaled K-gaps
  N1 <- sum(S_k > 0)
  N0 <- N_u - 1 - N1
  sum_qs <- sum(q_u * S_k)
  # Store the number of K-gaps, for use by nobs.kgaps()
  n_kgaps <- N0 + N1
  # Include censored inter-exceedance times?
  if (inc_cens) {
    n_kgaps <- n_kgaps + 2
    # censored inter-exceedance times and K-gaps
    T_u_cens <- c(exc_u[1] - 1, nx - exc_u[N_u])
    S_k_cens <- pmax(T_u_cens - k, 0)
    # N0, N1, sum of scaled K-gaps
    N1_cens <- sum(S_k_cens > 0)
    sum_s_cens <- sum(q_u * S_k_cens)
    # Add contributions.
    # Note: we divide N1_cens by two because a censored non-zero K-gap S_c
    # contributes theta exp(-theta q_u S_c) to the K-gaps likelihood,
    # whereas a non-censored non-zero K-gap contributes
    # theta^2 exp(-theta q_u S_c).
    # See equation (4.3) of Attalides (2015)
    N1 <- N1 + N1_cens / 2
    sum_qs <- sum_qs + sum_s_cens
  }
  return(list(N0 = N0, N1 = N1, sum_qs = sum_qs, n_kgaps = n_kgaps))
}

