# =================================== iwls ====================================
#
#' Iterated weighted least squares estimation of the extremal index
#'
#' Estimates the extremal index \eqn{\theta} using the iterated weighted least
#' squares method of Suveges (2007).  At the moment no estimates of
#' uncertainty are provided.
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u A numeric scalar.  Extreme value threshold applied to data.
#' @param maxit A numeric scalar.  The maximum number of iterations.
#' @details The iterated weighted least squares algorithm on page 46 of
#'   Suveges (2007) is used to estimate the value of the extremal index.
#'   This approach uses the time \emph{gaps} between successive exceedances
#'   in the data \code{data} of the threshold \code{u}.  The \eqn{i}th
#'   gap is defined as \eqn{T_i - 1}, where \eqn{T_i} is the difference in
#'   the occurrence time of exceedance \eqn{i} and exceedance \eqn{i + 1}.
#'   Therefore, threshold exceedances at adjacent time points produce a gap
#'   of zero.
#'
#'   The model underlying this approach is an exponential-point mas mixture
#'   for \emph{scaled gaps}, that is, gaps multiplied by the proportion of
#'   values in  \code{data} that exceed \code{u}.  Under this model
#'   scaled gaps are zero (`within-cluster' interexceedance times) with
#'   probability \eqn{1 - \theta} and otherwise (`between-cluster'
#'   interexceedance times) follow an exponential distribution with mean
#'   \eqn{1 / \theta}.
#'   The estimation method is based on fitting the `broken stick' model of
#'   Ferro (2003) to an exponential quantile-quantile plot of all of the
#'   scaled gaps.  Specifically, the broken stick is a horizontal line
#'   and a line with gradient \eqn{1 / \theta} which intersect at
#'   \eqn{(-\log\theta, 0)}{(-log \theta, 0)}.  The algorithm on page 46 of
#'   Suveges (2007) uses a weighted least squares minimization applied to
#'   the exponential
#'   part of this model to seek a compromise between the role of \eqn{\theta}
#'   as the proportion of interexceedance times that are between-cluster
#'   and the reciprocal of the mean of an exponential distribution for these
#'   interexceedance times.  The weights (see Ferro (2003)) are based on the
#'   variances of order statistics of a standard exponential sample: larger
#'   order statistics have larger sampling variabilities and therefore
#'   receive smaller weight than smaller order statistics.
#'
#'   Note that in step (1) of the algorithm on page 46 of Suveges there is a
#'   typo: \eqn{N_c + 1} should be \eqn{N}, where \eqn{N} is the number of
#'   threshold exceedances.  Also, the gaps are scaled as detailed above,
#'   not by their mean.
#' @return An object (a list) of class \code{"iwls", "exdex"} containing
#'     \item{\code{theta} }{The estimate of \eqn{\theta}.}
#'     \item{\code{conv} }{A convergence indicator: 0 indicates successful
#'       convergence; 1 indicates that \code{maxit} has been reached.}
#'     \item{\code{niter} }{The number of iterations performed.}
#'     \item{\code{n_gaps }}{The number of time gaps between successive
#'       exceedances.}
#'     \item{\code{call }}{The call to \code{iwls}.}
#' @references Suveges, M. (2007) Likelihood estimation of the extremal
#'   index. \emph{Extremes}, \strong{10}, 41-55.
#'   \url{https://doi.org/10.1007/s10687-007-0034-2}
#' @references Ferro, C.A.T. (2003) Statistical methods for clusters of
#'   extreme values. Ph.D. thesis, Lancaster University.
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{spm}} for estimation of the extremal index
#'   \eqn{\theta} using a semiparametric maxima method.
#' @examples
#' ### S&P 500 index
#'
#' u <- quantile(sp500, probs = 0.60)
#' theta <- iwls(sp500, u)
#' theta
#'
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = 0.90)
#' theta <- iwls(newlyn, u)
#' theta
#' @export
iwls <- function(data, u, maxit = 100) {
  Call <- match.call(expand.dots = TRUE)
  if (!is.numeric(u) || length(u) != 1) {
    stop("u must be a numeric scalar")
  }
  if (u >= max(data)) {
    stop("u must be less than max(data)")
  }
  # Calculate the quantities required to call iwls_fun()
  #
  # Sample size, positions, number and proportion of exceedances
  nx <- length(data)
  exc_u <- (1:nx)[data > u]
  N <- length(exc_u)
  # Inter-exceedances times, (largest first) 1-gaps, number of non-zero 1-gaps
  T_u <- diff(exc_u)
  S_1 <- pmax(T_u - 1, 0)
  n_gaps <- length(S_1)
  # Initial value of n_wls (the number of non-zero 1-gaps)
  n_wls <- length(S_1 > 0)
  # Sort the 1-gaps (largest to smallest) and scale by the sample proportion
  # of values that exceed u
  # [Bottom of page 45 of Suveges (2007), but not mentioned in the algorithm]
  S_1_sort <- sort(S_1, decreasing = TRUE)
  qhat <- N / nx
  S_1_sort <- S_1_sort * qhat
  # Standard exponential quantiles (based on ALL the N-1 1-gaps)
  exp_qs <- -log(1:(N - 1) / N)
  # Weights for the LS fit.  Larger value have large sampling variability and
  # therefore have smaller weights
  ws <-  rev(1 / cumsum(1 / (N:1) ^ 2))
  old_n_wls <- n_wls
  diff_n_wls <- 1
  niter <- 1L
  while (diff_n_wls != 0 & niter < maxit) {
    temp <- iwls_fun(n_wls, N, S_1_sort, exp_qs, ws, nx)
    n_wls <- temp$n_wls
    diff_n_wls <- n_wls - old_n_wls
    old_n_wls <- n_wls
    niter <- niter + 1L
  }
  conv <- ifelse(diff_n_wls > 0, 1, 0)
  n_wls <- temp$n_wls
  theta <- temp$theta
  res <- list(theta = theta, conv = conv, niter = niter, n_gaps = n_gaps,
              call = Call)
  class(res) <- c("iwls", "exdex")
  return(res)
}
