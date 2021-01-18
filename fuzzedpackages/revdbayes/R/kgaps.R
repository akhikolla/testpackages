# ================================= kgaps_post ================================

#' Random sampling from K-gaps posterior distribution
#'
#' Uses the \code{\link[rust]{rust}} package to simulate from the posterior
#' distribution of the extremal index \eqn{\theta} based on the K-gaps model
#' for threshold interexceedance times of Suveges and Davison (2010).
#'
#' @param n A numeric scalar. The size of posterior sample required.
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.  See Attalides (2015) for details.
#' @param alpha,beta Positive numeric scalars.  Parameters of a
#'   beta(\eqn{\alpha}, \eqn{\beta}) prior for \eqn{\theta}.
#' @param param A character scalar.  If \code{param = "logit"} (the default)
#'   then we simulate from the posterior distribution of
#'   \eqn{\phi = \log(\theta / (1-\theta))}{\phi = log(\theta / (1-\theta))}
#'   and then transform back to the
#'   \eqn{\theta}-scale.  If \code{param = "theta"} then we simulate
#'   directly from the posterior distribution of \eqn{\theta}, unless
#'   the sample K-gaps are all equal to zero or all positive, when we revert
#'   to \code{param = "logit"}.  This is to avoid sampling directly from a
#'   posterior with mode equal to 0 or 1.
#' @param use_rcpp A logical scalar.  If \code{TRUE} (the default) the
#'   rust function \code{\link[rust]{ru_rcpp}} is used for
#'   posterior simulation.  If \code{FALSE} the (slower) function
#'   \code{\link[rust]{ru}} is used.
#' @details A beta(\eqn{\alpha}, \eqn{\beta}) prior distribution is used for
#'   \eqn{\theta} so that the posterior from which values are simulated is
#'   proportional to
#'   \deqn{\theta ^ {2 N_1 + \alpha - 1} (1 - \theta) ^ {N_0 + \beta - 1}
#'     \exp\{- \theta q (S_0 + \cdots + S_N)\}.}{%
#'     \theta ^ (2 N_1 + \alpha - 1) * (1 - \theta) ^ (N_0 + \beta - 1) *
#'     exp(- \theta q (S_0 + \dots + S_N)).}
#'   See \code{\link{kgaps_stats}} for a description of the variables
#'   involved in the contribution of the likelihood to this expression.
#'
#'   The \code{\link[rust]{ru}} function in the \code{\link[rust]{rust}}
#'   package simulates from this posterior distribution using the
#'   generalised ratio-of-uniforms distribution.  To improve the probability
#'   of acceptance, and to ensure that the simulation will work even in
#'   extreme cases where the posterior density of \eqn{\theta} is unbounded as
#'   \eqn{\theta} approaches 0 or 1, we simulate from the posterior
#'   distribution of
#'   \eqn{\phi = \log(\theta / (1-\theta))}{\phi = log(\theta / (1-\theta))}
#'   and then transform back to the \eqn{\theta}-scale.
#' @return An object (list) of class \code{"evpost"}, which has the same
#'   structure as an object of class \code{"ru"} returned from
#'   \code{\link[rust]{ru}}.
#'   In addition this list contains
#'   \itemize{
#'     \item{\code{model}:} The character scalar \code{"kgaps"}.
#'     \item{\code{thresh}:} The argument \code{thresh}.
#'     \item{\code{ss}:} The sufficient statistics for the K-gaps likelihood,
#'       as calculated by \code{\link{kgaps_stats}}.
#'   }
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#' @seealso \code{\link{kgaps_mle}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @seealso \code{\link{kgaps_stats}} for the calculation of sufficient
#'   statistics for the K-gaps model.
#' @seealso \code{\link[rust]{ru}} for the form of the object returned by
#'   \code{kgaps_post}.
#' @examples
#' thresh <- quantile(newlyn, probs = 0.90)
#' k_postsim <- kgaps_post(newlyn, thresh)
#' plot(k_postsim)
#' @export
kgaps_post <- function(data, thresh, k = 1, n = 1000, inc_cens = FALSE,
                       alpha = 1, beta = 1, param = c("logit", "theta"),
                       use_rcpp = TRUE) {
  if (!is.numeric(thresh) || length(thresh) != 1) {
    stop("thresh must be a numeric scalar")
  }
  if (thresh >= max(data)) {
    stop("thresh must be less than max(data)")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("k must be a numeric scalar")
  }
  if (k < 1) {
    stop("k must be no smaller than 1.")
  }
  if (alpha <= 0 | beta <= 0) {
    stop("alpha and beta must be positive.")
  }
  param <- match.arg(param)
  # Calculate the MLE and get the sufficient statistics
  mle_list <- kgaps_mle(data, thresh, k, inc_cens)
  theta_mle <- mle_list$theta_mle
  ss <- mle_list$ss
  # Set an initial value for theta, and perhaps phi. We do this by noting that
  # (a) the Jacobian of the transformation from theta to phi multiplies the
  # likelihood by a factor of theta * (1 - theta), and (b) the beta prior
  # includes a factor of theta ^ (alpha - 1) * (1 - theta) ^ (beta - 1).
  # Therefore, we can find the MAP of theta by solving a quadratic equation
  # for theta, in the same way that we do to find the MLE for theta.
  #
  theta_init <- kgaps_quad_solve(ss$N0 + alpha, ss$N1 + beta, ss$sum_qs)
  # Set essential arguments to ru()
  # ... including the K-gaps posterior distribution for theta
  if (use_rcpp) {
    post_ptr <- kgaps_logpost_xptr("kgaps")
    for_post <- c(ss, list(alpha = alpha, beta = beta))
    for_ru <- list(logf = post_ptr, pars = for_post, n = n)
  } else {
    logpost <- function(theta, ss) {
      loglik <- do.call(kgaps_loglik, c(list(theta = theta), ss))
      if (is.infinite(loglik)) return(loglik)
      # Add beta(alpha, beta) prior
      logprior <- stats::dbeta(theta, alpha, beta, log = TRUE)
      return(loglik + logprior)
    }
    for_ru <- list(logf = logpost, ss = ss, n = n)
  }
  # Sample on the logit scale phi = log(theta / (1 - theta)) ?
  if (ss$N0 == 0 || ss$N1 == 0) {
    param = "logit"
  }
  if (param == "logit") {
    # Transformation, Jacobian and initial estimate
    if (use_rcpp) {
      phi_to_theta <- phi_to_theta_xptr("kgaps")
      log_j <- log_j_xptr("kgaps")
    } else {
      phi_to_theta <- function(phi) {
        ephi <- exp(phi)
        return(ephi / (1 + ephi))
      }
      log_j <- function(theta) {
        return(-log(theta) - log(1 - theta))
      }
    }
    phi_init <- log(theta_init / (1 - theta_init))
    trans_list <- list(phi_to_theta = phi_to_theta, log_j = log_j)
    for_ru <- c(for_ru, list(init = phi_init, trans = "user"), trans_list)
  } else {
    for_ru <- c(for_ru, list(init = theta_init, trans = "none"))
  }
  ru_fn <- ifelse(use_rcpp, rust::ru_rcpp, rust::ru)
  temp <- do.call(ru_fn, for_ru)
  temp$model <- "kgaps"
  temp$thresh <- thresh
  temp$ss <- ss
  class(temp) <- "evpost"
  return(temp)
}

# ================================= kgaps_mle =================================
#
#' Maximum likelihood estimation for the K-gaps model
#'
#' Calculates maximum likelihood estimates of the extremal index \eqn{\theta}
#' based on the K-gaps model for threshold inter-exceedances times of
#' Suveges and Davison (2010).
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.  See Attalides (2015) for details.
#' @param conf  A numeric scalar.  If \code{conf} is supplied then a
#'   \code{conf}\% likelihood-based confidence interval for \eqn{\theta} is
#'   estimated.
#' @details The maximum likelihood estimate of the extremal index \eqn{\theta}
#'   under the K-gaps model of Suveges and Davison (2010) is calculated.
#'   If \code{inc_cens = TRUE} then information from censored inter-exceedance
#'   times is included in the likelihood to be maximised, following
#'   Attalides (2015).  The form of the log-likelihood is given in the
#'   \strong{Details} section of \code{\link{kgaps_stats}}.
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#' @return A list containing
#'   \itemize{
#'     \item {\code{theta_mle} : } {The maximum likelihood estimate (MLE) of
#'       \eqn{\theta}.}
#'     \item {\code{theta_se} : } {The estimated standard error of the MLE.}
#'     \item {\code{theta_ci} : } {(If \code{conf} is supplied) a numeric
#'       vector of length two giving lower and upper confidence limits for
#'       \eqn{\theta}.}
#'     \item {\code{ss} : } {The list of summary statistics returned from
#'       \code{\link{kgaps_stats}}.}
#'   }
#' @seealso \code{\link{kgaps_stats}} for the calculation of sufficient
#'   statistics for the K-gaps model.
#' @examples
#' thresh <- quantile(newlyn, probs = 0.90)
#' # MLE and SE only
#' kgaps_mle(newlyn, thresh)
#' # MLE, SE and 95% confidence interval
#' kgaps_mle(newlyn, thresh, conf = 95)
#' @export
kgaps_mle <- function(data, thresh, k = 1, inc_cens = FALSE, conf = NULL) {
  if (!is.numeric(thresh) || length(thresh) != 1) {
    stop("thresh must be a numeric scalar")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("k must be a numeric scalar")
  }
  if (thresh >= max(data)) {
    stop("thresh must be less than max(data)")
  }
  # Calculate sufficient statistics
  ss <- kgaps_stats(data, thresh, k, inc_cens)
  # If N0 = 0 then all exceedances occur singly (all K-gaps are positive)
  # and the likelihood is maximised at theta = 1.
  N0 <- ss$N0
  # If N1 = 0 then we are in the degenerate case where there is one cluster
  # (all K-gaps are zero) and the likelihood is maximised at theta = 0.
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
  if (is.null(conf)) {
    return(list(theta_mle = theta_mle, theta_se = theta_se, ss = ss))
  }
  conf_int <- kgaps_conf_int(theta_mle, ss, conf = conf)
  return(list(theta_mle = theta_mle, theta_se = theta_se, theta_ci = conf_int,
              ss = ss))
}

# ================================ kgaps_stats ================================

#' Sufficient statistics for the K-gaps model
#'
#' Calculates sufficient statistics for the K-gaps model for the extremal index
#' \eqn{\theta}.
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param thresh A numeric scalar.  Extreme value threshold applied to data.
#' @param k A numeric scalar.  Run parameter \eqn{K}, as defined in Suveges and
#'   Davison (2010).  Threshold inter-exceedances times that are not larger
#'   than \code{k} units are assigned to the same cluster, resulting in a
#'   \eqn{K}-gap equal to zero.  Specifically, the \eqn{K}-gap \eqn{S}
#'   corresponding to an inter-exceedance time of \eqn{T} is given by
#'   \eqn{S = max(T - K, 0)}.
#' @param inc_cens A logical scalar indicating whether or not to include
#'   contributions from censored inter-exceedance times relating to the
#'   first and last observation.  See Attalides (2015) for details.
#' @details The sample K-gaps are
#'   \eqn{S_0, S_1, \ldots, S_{N-1}, S_N}{S_0, S_1, ..., S_(N-1), S_N}, where
#'   \eqn{S_1, \ldots, S_{N-1}}{S_1, ..., S_(N-1)} are uncensored and \eqn{S_0}
#'   and \eqn{S_N} are censored.  Under the assumption that the K-gaps are
#'   independent, the log-likelihood of the K-gaps model is given by
#'   \deqn{l(\theta; S_0, \ldots, S_N) = N_0 \log(1 - \theta) +
#'     2 N_1 \log \theta - \theta q (S_0 + \cdots + S_N),}{%
#'     l(\theta; S_0, ..., S_N) = N_0 log(1 - \theta) + 2 N_1 log \theta -
#'     \theta q (S_0 + ... + S_N),}
#'    where \eqn{q} is the threshold exceedance probability,
#'    \eqn{N_0} is the number of sample K-gaps that are equal to zero and
#'    (apart from an adjustment for the contributions of \eqn{S_0} and
#'    \eqn{S_N}) \eqn{N_1} is the number of positive sample K-gaps.
#'    Specifically, \eqn{N_1} is equal to the number of
#'    \eqn{S_1, \ldots, S_{N-1}}{S_1, ..., S_(N-1)}
#'    that are positive plus \eqn{(I_0 + I_N) / 2}, where \eqn{I_0 = 1} if
#'    \eqn{S_0} is greater than zero and similarly for \eqn{I_N}.
#'    The differing treatment of uncensored and censored K-gaps reflects
#'    differing contributions to the likelihood.
#'    For full details see Suveges and Davison (2010) and Attalides (2015).
#' @return A list containing the sufficient statistics, with components
#'   \itemize{
#'     \item {\code{N0} : } {the number of zero K-gaps}
#'     \item {\code{N1} : } {contribution from non-zero K-gaps (see
#'       \strong{Details})}
#'     \item {\code{sum_qs} : } {the sum of the (scaled) K-gaps, i.e.
#'       \eqn{q (S_0 + \cdots + S_N)}{q (S_0 + ... + S_N)}, where \eqn{q} is
#'       estimated by the proportion of threshold exceedances.}
#'   }
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @references Attalides, N. (2015) Threshold-based extreme value modelling,
#'   PhD thesis, University College London.
#' @seealso \code{\link{kgaps_mle}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the K-gaps model.
#' @examples
#' u <- quantile(newlyn, probs = 0.90)
#' kgaps_stats(newlyn, u)
#' @export
kgaps_stats <- function(data, thresh, k = 1, inc_cens = FALSE) {
  if (any(is.na(data))) {
    stop("No missing values are allowed in ''data''")
  }
  if (!is.numeric(thresh) || length(thresh) != 1) {
    stop("thresh must be a numeric scalar")
  }
  if (thresh >= max(data)) {
    stop("thresh must be less than max(data)")
  }
  if (!is.numeric(k) || length(k) != 1) {
    stop("k must be a numeric scalar")
  }
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > thresh]
  N_u <- length(exc_u)
  q_u <- N_u / nx
  # Inter-exceedances times and K-gaps
  T_u <- diff(exc_u)
  S_k <- pmax(T_u - k,0)
  # N0, N1, sum of scaled K-gaps
  N1 <- sum(S_k > 0)
  N0 <- N_u - 1 - N1
  sum_qs <- sum(q_u * S_k)
  # Include censored inter-exceedance times?
  if (inc_cens) {
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
  return(list(N0 = N0, N1 = N1, sum_qs = sum_qs))
}

# =============================== kgaps_loglik ================================

kgaps_loglik <- function(theta, N0, N1, sum_qs){
  if (theta < 0 || theta > 1) {
    return(-Inf)
  }
  loglik <- 0
  if (N1 > 0) {
    loglik <- loglik + 2 * N1 * log(theta) - sum_qs * theta
  }
  if (N0 > 0) {
    loglik <- loglik + N0 * log(1 - theta)
  }
  return(loglik)
}

# ============================== kgaps_conf_int ===============================

kgaps_conf_int <- function(theta_mle, ss, conf = 95) {
  cutoff <- stats::qchisq(conf / 100, df = 1)
  theta_list <- c(list(theta = theta_mle), ss)
  max_loglik <- do.call(kgaps_loglik, theta_list)
  ob_fn <- function(theta) {
    theta_list$theta <- theta
    loglik <- do.call(kgaps_loglik, theta_list)
    return(2 * (max_loglik - loglik) - cutoff)
  }
  ci_low <- 0
  ci_up <- 1
  if (ss$N1 > 0) {
    ci_low <- stats::uniroot(ob_fn, c(0, theta_mle))$root
  }
  if (ss$N0 > 0) {
    ci_up <- stats::uniroot(ob_fn, c(theta_mle, 1))$root
  }
  return(c(ci_low, ci_up))
}

# ============================== kgaps_quad_solve =============================

kgaps_quad_solve <- function(N0, N1, sum_qs) {
  aa <- sum_qs
  bb <- -(N0 + 2 * N1 + sum_qs)
  cc <- 2 * N1
  qq <- -(bb - sqrt(bb ^ 2 - 4 * aa * cc)) / 2
  theta_mle <- cc / qq
  return(theta_mle)
}

