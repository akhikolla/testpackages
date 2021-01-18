# =============================== spm() ===================================== #

#' Semiparametric maxima estimator of the extremal index
#'
#' Estimates the extremal index \eqn{\theta} using a semiparametric block
#' maxima estimator of Northrop (2015) (\code{"N2015"}) and a variant of this
#' estimator studied by Berghaus and Bucher (2018) (\code{"BB2018"}), using
#' both sliding (overlapping) block maxima and disjoint (non-overlapping) block
#' maxima.  A simple modification (subtraction of \eqn{1 / b}, where \eqn{b} is
#' the block size) of the Berghaus and Bucher (2018) estimator
#' (\code{"BB2018b"}) is also calculated.  Estimates of uncertainty are
#' calculated using the asymptotic theory developed by Berghaus and Bucher
#' (2018).
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param b A numeric scalar.  The block size.
#' @param bias_adjust A character scalar.  Is bias-adjustment of the
#'   raw estimate of \eqn{\theta} performed using the bias-reduced
#'   estimator (\code{bias_adjust = "BB3"}), derived in Section 5 of
#'   Berghaus and Bucher (2018); or a simpler version
#'   (\code{bias_adjust = "BB1"}), in which the raw estimate is multiplied by
#'   \eqn{(k-1) / k}, where \eqn{k} is the number of blocks; or the
#'   bias-adjustment of the empirical distribution function used to calculate
#'   the estimate, as detailed in Section 2 of Northrop (2015).  When disjoint
#'   maxima are used \code{bias_adjust = "BB1"} and \code{bias_adjust = "N"}
#'   give identical estimates of the Berghaus and Bucher (2018) variant,
#'   as explained at the end of Section 5 of Berghaus and Bucher (2018).
#'   If \code{bias_adjust = "none"} then no bias-adjustment is performed.
#' @param constrain A logical scalar.  If \code{constrain = TRUE} then
#'   any estimates that are greater than 1 are set to 1,
#'   that is, they are constrained to lie in (0, 1].  This is carried out
#'   \emph{after} any bias-adjustment.  Otherwise,
#'   estimates that are greater than 1 may be obtained.
#' @param varN A logical scalar.  If \code{varN = TRUE} then the estimation
#'   of the sampling variance of the Northrop (2015) estimator is tailored
#'   to that estimator.  Otherwise, the sampling variance derived in
#'   Berghaus and Bucher (2018) is used.
#'   See \strong{Details} for further information.
#' @param which_dj A character scalar.  Determines which set of disjoint
#'   maxima are used to calculate an estimate of \eqn{\theta}: \code{"first"},
#'   only the set whose first block starts on the first observation in
#'   \code{x}; \code{"last"}, only the set whose last block ends on the last
#'   observation in \code{x}.
#' @details The extremal index \eqn{\theta} is estimated using the
#'   semiparametric maxima estimator of Northrop (2015) and variant
#'   of this studied by Berghaus and Bucher (2018).  In each case a sample
#'   of 'data' is derived from the input data \code{data}, based on the
#'   empirical distribution function of these data evaluated at the
#'   maximum values of  of blocks of \code{b} contiguous values in \code{data}.
#'
#'   The estimators are based on an assumption that these 'data' are sampled
#'   approximately from an exponential distribution with mean \eqn{1/\theta}.
#'   For details see page 2309 of Berghaus and Bucher (2018), where the
#'   'data' for the Northrop (2015) estimator are denoted \eqn{Y} and
#'   those for the Berghaus and Bucher (2018) are denoted \eqn{Z}.
#'   For convenience, we will refer to these values as the
#'   \eqn{Y}-data and the \eqn{Z}-data.
#'
#'   The approximate nature of the model for the \eqn{Y}-data arises from the
#'   estimation of the distribution function \eqn{F}.  A further
#'   approximation motivates the use of the \eqn{Z}-data.  If \eqn{F} is known
#'   then the variable \eqn{Z / b} has a beta(1, \eqn{b\theta}) distribution,
#'   so that that is, \eqn{Z} has mean \eqn{1 / (\theta + 1/b)}.
#'   Therefore, an exponential distribution with mean \eqn{1 / (\theta + 1/b)}
#'   may provide a better approximate model, which provides the motivation for
#'   subtracting \eqn{1 / b} from the Berghaus and Bucher (2018) estimator.
#'   Indeed, the resulting estimates are typically close to those
#'   of the Northrop (2015) estimator.
#'
#'   If \code{sliding = TRUE} then the function uses sliding block maxima,
#'   that is, the largest value observed in \emph{all}
#'   (\code{length(data) - b + 1}) blocks of \code{b} observations.
#'   If \code{sliding = FALSE} then disjoint block maxima, that is, the
#'   largest values in (\code{floor(n / b)}) disjoint blocks of \code{b}
#'   observations, are used.
#'
#'   Estimation of the sampling variances of the estimators is based
#'   on Proposition 4.1 on page 2319 of Berghaus and Bucher (2018).
#'   For the Northrop (2015) variant the user has the choice either to
#'   use the sampling variance based on the Berghaus and Bucher (2018)
#'   estimator, i.e. the \eqn{Z}-data (\code{varN = FALSE}) or an analogous
#'   version tailored to the Northrop (2015) estimator that uses the
#'   \eqn{Y}-data (\code{varN = TRUE}).
#'
#'   The estimates of the sampling variances of the sliding blocks estimators
#'   are inferred from those of the disjoint blocks estimators (see page 2319
#'   of Berhaus and Bucher (2018)).  The calculation of the latter uses a set
#'   of disjoint block maxima.  If \code{length(data)} is not an integer
#'   multiple of \code{b} then there will be more than one set of these, and
#'   all are equally valid.  In this event we perform the calculation for all
#'   such sets and use the mean of the resulting estimates.  This reduces the
#'   sampling variability of the estimates at the expense of slowing down
#'   the calculation somewhat, particularly if \code{b} is small.  This may
#'   become apparent when calling \code{spm} repeatedly in
#'   \code{\link{choose_b}}.
#'
#'   This estimator of the sampling variance of the sliding blocks estimator is
#'   not constrained to be positive: a negative estimate may result if the
#'   block size is small.  In this event
#'   \strong{no warning will be given until the returned object is printed} and,
#'   for the affected estimator (\code{"N2015"} or \code{"BB2018/BB2018b"}),
#'     \itemize{
#'       \item{the corresponding estimated standard errors using sliding
#'             blocks will be missing in \code{se_sl} in the returned object,}
#'       \item{if \code{bias_adjust == "BB3"} then bias-adjustment
#'             based on \code{bias_adjust == "BB1"} will instead be performed
#'             when using sliding blocks, because the former relies on the
#'             estimated variances of the estimators.}
#'     }
#'  Similarly, bias adjustment under \code{adjust = "BB3"} and/or subtraction
#'  of \eqn{1 / b} in the \code{"BB2018b"} case may, rare cases, produce a
#'  negative estimate of \eqn{\theta}.  In these instances an estimate of
#'  zero is returned, but the values returned in \code{bias_dj} and
#'  \code{bias_sl} are not changed.
#' @return A list of class \code{c("spm", "exdex")} containing the
#'   components listed below.  The components that are vectors are
#'   labelled to indicate the estimator to which the constituent values
#'   relate: \code{"N2015"} for Northrop (2015), \code{"BB2018"} for
#'   Berghaus and Bucher (2018) and \code{"BB2018b"} for the modified version.
#'   \item{theta_sl, theta_dj}{ Vectors containing the estimates of
#'     \eqn{\theta} resulting from sliding maxima and disjoint maxima
#'     respectively.}
#'   \item{se_sl, se_dj}{The estimated standard errors associated
#'     with the estimates in \code{theta_sl} and \code{theta_dj}.
#'     The values for \code{"BB2018"} and \code{"BB2018b"} are identical.}
#'   \item{bias_sl, bias_dj}{The respective values of the
#'       bias-adjustment applied to the raw estimates, that is, the values
#'       subtracted from the raw estimates.  For estimator
#'       \code{BB2018b} this includes a contribution for the subtraction
#'       of \code{1 / b}.
#'       If \code{bias_adjust = "N"} or \code{"none"} then
#'       \code{bias_sl} and \code{bias_dj} are \code{c(0, 0, 1 / b)}.}
#'   \item{raw_theta_sl, raw_theta_dj}{ Vectors containing the raw estimates
#'     of \eqn{\theta}, prior to any bias_adjustment.}
#'   \item{uncon_theta_sl, uncon_theta_dj}{The (bias_adjusted) estimates of
#'     \eqn{\theta} before the constraint that they lie in (0, 1] has been
#'     applied.}
#'   \item{data_sl, data_dj}{Matrices containing the \eqn{Y}-data and
#'     \eqn{Z}-data for the sliding an disjoint maxima respectively.
#'     The first columns are the \eqn{Y}-data, the second columns the
#'     \eqn{Z}-data.}
#'   \item{sigma2dj, sigma2dj_for_sl}{Estimates of the variance
#'    \eqn{\sigma_{{\rm dj}}^2}{\sigma^2_dj}
#'    defined on pages 2314-2315 of Berghaus and Bucher (2018).
#'    The form of the estimates is given on page 2319.
#'    \code{sigma2dj} is used in estimating the standard error \code{se_dj},
#'    \code{sigma2dj_for_sl} in estimating the standard error \code{se_sl}.}
#'   \item{sigma2sl}{Estimates of the variance
#'     \eqn{\sigma_{{\rm sl}}^2}{\sigma^2_sl}.
#'     defined on pages 2314-2315 of Berghaus and Bucher (2018).
#'     The form of the estimates is given on page 2319.
#'     \code{sigma2dj_for_sl} is used to estimate
#'     \eqn{\sigma_{{\rm dj}}^2}{\sigma^2_dj} for this purpose.}
#'   \item{b}{The input value of \code{b}.}
#'   \item{bias_adjust}{The input value of \code{bias_adjust}.}
#'   \item{call}{The call to \code{spm}.}
#' @seealso \code{\link{confint.spm}} to estimate confidence intervals
#'   for \eqn{theta}.
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{iwls}}: iterated weighted least squares estimator.
#' @references Northrop, P. J. (2015) An efficient semiparametric maxima
#' estimator of the extremal index. \emph{Extremes} \strong{18}(4), 585-603.
#' \url{https://doi.org/10.1007/s10687-015-0221-5}
#' @references Berghaus, B., Bucher, A. (2018) Weak convergence of a pseudo
#' maximum likelihood estimator for the extremal index. \emph{Ann. Statist.}
#' \strong{46}(5), 2307-2335. \url{https://doi.org/10.1214/17-AOS1621}
#' @examples
#' ### Newlyn sea surges
#' theta <- spm(newlyn, 20)
#' theta
#' summary(theta)
#'
#' ### S&P 500 index
#' theta <- spm(sp500, 100)
#' theta
#' summary(theta)
#' @export
spm <- function(data, b, bias_adjust = c("BB3", "BB1", "N", "none"),
                constrain = TRUE, varN = TRUE,
                which_dj = c("last", "first")) {
  Call <- match.call(expand.dots = TRUE)
  #
  # Check inputs
  #
  if (missing(data) || length(data) == 0L || mode(data) != "numeric") {
    stop("'data' must be a non-empty numeric vector")
  }
  if (any(!is.finite(data))) {
    stop("'data' contains missing or infinite values")
  }
  if (is.matrix(data)) {
    stop("'data' must be a vector")
  }
  data <- as.vector(data)
  if (!is.numeric(b) || length(b) != 1) {
    stop("'b' must be a numeric scalar (specifically, a positive integer)")
  }
  if (b < 1) {
    stop("'b' must be no smaller than 1")
  }
  bias_adjust <- match.arg(bias_adjust)
  if (!is.logical(constrain) || length(constrain) != 1) {
    stop("'constrain' must be a logical scalar")
  }
  if (!is.logical(varN) || length(varN) != 1) {
    stop("'varN' must be a logical scalar")
  }
  which_dj <- match.arg(which_dj)
  # Find the number of (disjoint) blocks
  k_n <- floor(length(data) / b)
  if (k_n < 1) {
    stop("b is too large: it is larger than length(data)")
  }
  #
  # Estimate sigma2_dj based on Section 4 of Berghaus and Bucher (2018)
  # We require the disjoint maxima to do this.  If sliding = TRUE then
  # pass these to spm_sigmahat_dj using the dj_maxima argument
  # Only do this is b_ok = TRUE.
  # Otherwise, just calculate point estimates of theta
  # At this point these estimates have not been bias-adjusted, unless
  # bias_adjust = "N".
  #
  # Find all sets of maxima of disjoint blocks of length b
  all_max <- all_max_rcpp(data, b, which_dj = "all")
  # Here, m is the number of xs that contribute to each set of disjoint block
  # maxima
  m <- nrow(all_max$xd)
  res <- cpp_sigma2hat_dj(all_max = all_max, b = b, kn = k_n, m = m,
                          bias_adjust = bias_adjust, which_dj = which_dj)
  est_names <- c("N2015", "BB2018")
  # In res theta_dj ... are 2x1 matrices.  Convert them to named vectors.
  res$theta_dj <- as.vector(res$theta_dj)
  names(res$theta_dj) <- est_names
  res$sigma2dj <- as.vector(res$sigma2dj)
  names(res$sigma2dj) <- est_names
  res$sigma2dj_for_sl <- as.vector(res$sigma2dj_for_sl)
  names(res$sigma2dj_for_sl) <- est_names
  colnames(res$data_dj) <- est_names
  # Sliding maxima
  Fhaty <- ecdf2(all_max$xs, all_max$ys)
  # Avoid over-writing the `disjoint' sample size k_n: it is needed later
  k_n_sl <- length(all_max$ys)
  # Now, m is the number of xs that contribute to the sliding maxima
  m <- length(all_max$xs)
  const <- -log(m - b + k_n_sl)
  if (bias_adjust == "N") {
    Fhaty <- (m * Fhaty - b) / (m - b)
  }
  res$theta_sl <- c(-1 / mean(b * log0const(Fhaty, const)),
                    1 / (b * mean(1 - Fhaty)))
  names(res$theta_sl) <- est_names
  #
  # Add the values of the Y-data and the Z-data to the output
  res$data_sl <- cbind(N2015 = -b * log(Fhaty), BB2018 = b * (1 - Fhaty))
  #
  # Estimate the sampling variances of the estimators
  #
  res$sigma2sl <- res$sigma2dj_for_sl - (3 - 4 * log(2)) / res$theta_sl ^ 2
  # res$sigma2sl could contain non-positive values
  # If it does then replace them with NA
  res$sigma2sl[res$sigma2sl <= 0] <- NA
  indexN <- ifelse(varN, 2, 1)
  if (varN) {
    index <- 1:2
  } else {
    index <- c(2, 2)
  }
  res$se_dj <- res$theta_dj ^ 2 * sqrt(res$sigma2dj[index] / k_n)
  res$se_sl <- res$theta_sl ^ 2 * sqrt(res$sigma2sl[index] / k_n)
  #
  # Store the raw (not bias-adjusted) estimates
  #
  res$raw_theta_dj <- res$theta_dj
  res$raw_theta_sl <- res$theta_sl
  #
  # Perform BB2018 bias-adjustment if required
  #
  if (bias_adjust == "BB3") {
    res$bias_dj <- res$theta_dj / k_n + res$theta_dj ^ 3 * res$sigma2dj / k_n
    res$theta_dj <- res$theta_dj - res$bias_dj
    BB3adj_sl <- res$theta_sl / k_n + res$theta_sl ^ 3 * res$sigma2sl / k_n
    BB1adj_sl <- res$theta_sl / k_n
    res$bias_sl <- ifelse(is.na(res$se_sl), BB1adj_sl, BB3adj_sl)
    res$theta_sl <- res$theta_sl - res$bias_sl
  } else if (bias_adjust == "BB1") {
    res$bias_dj <- res$theta_dj / k_n
    res$theta_dj <- res$theta_dj - res$bias_dj
    res$bias_sl <- res$theta_sl / k_n
    res$theta_sl <- res$theta_sl - res$bias_sl
  } else {
    res$bias_dj <- res$bias_sl <- c(N2015 = 0, BB2018 = 0)
  }
  #
  # Add estimates, bias and standard errors for an estimator that results from
  # applying a further bias-adjustment to the BB2018 estimator.  The adjustment
  # is to subtract 1/b.  This stems from the fact that the expected value of
  # the BB2018 variable Z is 1/(theta + 1/b), not 1/theta.
  #
  # Estimates
  res$theta_dj <- c(res$theta_dj, res$theta_dj["BB2018"] - 1 / b)
  names(res$theta_dj)[3] <- "BB2018b"
  res$theta_sl <- c(res$theta_sl, res$theta_sl["BB2018"] - 1 / b)
  names(res$theta_sl)[3] <- "BB2018b"
  # Bias adjustment
  if (bias_adjust == "BB3" || bias_adjust == "BB1") {
    res$bias_dj <- c(res$bias_dj, res$bias_dj["BB2018"] + 1 / b)
    res$bias_sl <- c(res$bias_sl, res$bias_sl["BB2018"] + 1 / b)
  } else {
    res$bias_dj <- c(res$bias_dj, 1 / b)
    res$bias_sl <- c(res$bias_sl, 1 / b)
  }
  names(res$bias_dj)[3] <- "BB2018b"
  names(res$bias_sl)[3] <- "BB2018b"
  # Standard errors
  res$se_dj <- c(res$se_dj, res$se_dj["BB2018"])
  names(res$se_dj)[3] <- "BB2018b"
  res$se_sl <- c(res$se_sl, res$se_sl["BB2018"])
  names(res$se_sl)[3] <- "BB2018b"
  # Save the unconstrained estimates, so that they can be returned
  res$uncon_theta_dj <- res$theta_dj
  res$uncon_theta_sl <- res$theta_sl
  # Constrain to (0, 1] if required
  if (constrain) {
    res$theta_dj <- pmin(res$theta_dj, 1L)
    res$theta_sl <- pmin(res$theta_sl, 1L)
  }
  # Constrain estimates to be non-negative. BB3 bias-adjustment could result in
  # negative estimates as could subtracting 1/b from BB2018 to form BB2018b.
  res$theta_dj <- pmax(res$theta_dj, 0L)
  res$theta_sl <- pmax(res$theta_sl, 0L)
  #
  res$bias_adjust <- bias_adjust
  res$b <- b
  res$call <- Call
  #
  class(res) <- c("spm", "exdex")
  return(res)
}
