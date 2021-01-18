# ================================= kgaps_imt =================================
#
#' Information matrix test under the \eqn{K}-gaps model
#'
#' Performs the information matrix test (IMT) of Suveges and Davison (2010) to
#' diagnose misspecification of the \eqn{K}-gaps model
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u,k Numeric vectors.  \code{u} is a vector of
#'   extreme value thresholds applied to data.  \code{k} is a vector of values
#'   of the run parameter \eqn{K}, as defined in Suveges and Davison (2010).
#'   See \code{\link{kgaps}} for more details.
#' @details The IMT is performed a over grid of all
#'   combinations of threshold and \eqn{K} in the vectors \code{u}
#'   and \code{k}.  If the estimate of \eqn{\theta} is 0 then the
#'   IMT statistic, and its associated \eqn{p}-value will be \code{NA}.
#'
#'   For details of the IMT see Suveges and Davison
#'   (2010).  There are some typing errors on pages 18-19 that have been
#'   corrected in producing the code: the penultimate term inside \code{{...}}
#'   in the middle equation on page 18 should be \eqn{(c_j(K))^2}, as should
#'   the penultimate term in the first equation on page 19; the \code{{...}}
#'   bracket should be squared in the 4th equation on page 19; the factor
#'   \eqn{n} should be \eqn{N-1} in the final equation on page 19.
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @return An object (a list) of class \code{c("kgaps_imt", "exdex")}
#'   containing
#'   \item{imt }{A \code{length(u)} by \code{length(k)} numeric matrix.
#'     Column i contains, for K = \code{k[i]}, the values of the
#'     information matrix test statistic for the set of thresholds in
#'     \code{u}.  The column names are the values in code{k}.
#'     The row names are the approximate empirical percentage quantile levels
#'     of the thresholds in \code{u}.}
#'   \item{p }{A \code{length(u)} by \code{length(k)} numeric matrix
#'     containing the corresponding \eqn{p}-values for the test.}
#'   \item{theta }{A \code{length(u)} by \code{length(k)} numeric matrix
#'     containing the corresponding estimates of \eqn{\theta}.}
#'   \item{u,k }{The input \code{u} and \code{k}.}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @examples
#' u <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
#' imt <- kgaps_imt(newlyn, u, k = 1:5)
#' @export
kgaps_imt <- function(data, u, k = 1) {
  # Function to return only the MLE of theta
  mle_only <- function(k, data, u) {
    return(kgaps(data, u, k, inc_cens = FALSE)$theta)
  }
  theta <- T_mat <- p_mat <- NULL
  n_u <- length(u)
  n_k <- length(k)
  # Beginning of loop over all thresholds ----------
  for (iu in 1:n_u) {
    the_u <- u[iu]
    # Calculate the MLE of theta for each value of k
    thetahat <- vapply(k, mle_only, 0, data = data, u = the_u)
    # sample size of x
    nx <- length(data)
    # positions of exceedances of u
    exc_u <- (1:nx)[data > the_u]
    # number of exceedances
    n_u <- length(exc_u)
    # proportion of values that exceed u
    q_u <- n_u / nx
    # inter-exceedance times
    T_u <- diff(exc_u)
    #
    # Create a length(T) by length(k) matrix with column i containing
    # the values of S_k for k=k[i] ...
    #
    # One column for each value of k
    S_k <- sapply(k, function(k) pmax(T_u - k, 0))
    c_mat <- q_u * S_k
    theta_mat <- matrix(thetahat, ncol = n_k, nrow = nrow(S_k), byrow = TRUE)
    ld <- ifelse(c_mat == 0, -1 / (1 - theta_mat), 2 / theta_mat) - c_mat
    neg_ldd <- ifelse(c_mat == 0, 1 / (1 - theta_mat) ^ 2, 2 / theta_mat ^ 2)
    In <- colMeans(neg_ldd)
    Jn <- colMeans(ld ^ 2)
    Dn <- Jn - In
    dc <- ld ^ 2 - neg_ldd
    dcd <- 4 * c_mat / theta_mat ^ 2 + ifelse(c_mat == 0, 0, -4 / theta_mat ^ 3)
    # Force NA, rather than NA, in cases where thetahat = 0
    dcd[is.nan(dcd)] <- NA
    Dnd <- colMeans(dcd)
    # Multiply the columns of ld by the corresponding elements of Dnd / In
    temp <- ld * rep(Dnd / In, rep(nrow(ld), ncol(ld)))
    Vn <- colMeans((dc - temp) ^ 2)
    test_stats <- (n_u - 1) * Dn ^ 2 / Vn
    pvals <- stats::pchisq(test_stats, df = 1, lower.tail = FALSE)
    theta <- rbind(theta, thetahat)
    T_mat <- rbind(T_mat, test_stats)
    p_mat <- rbind(p_mat, pvals)
  }
  # End of loop over thresholds ----------
  colnames(T_mat) <- colnames(p_mat) <- colnames(theta) <- k
  if (is.null(names(u))) {
    u_ps <- round(100 * sapply(u, function(x) mean(data < x)))
  } else {
    u_ps <- as.numeric(substr(names(u), 1, nchar(names(u),
                                                     type = "c") - 1))
  }
  rownames(T_mat) <- rownames(p_mat) <- rownames(theta) <- u_ps
  res <- list(imt = T_mat, p = p_mat, theta = theta, u = u, k = k)
  class(res) <- c("kgaps_imt", "exdex")
  return(res)
}
