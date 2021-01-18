#' Estimate quantiles of distribution
#'
#' @param p A vector of probabilities.
#' @param x The x-axis values corresponding to the cumulative distribution function.
#' @param cdf A cumulative distributions function, i.e. output of \code{\link{estCdf}}.
#'
#' @details Quantiles are obtained in the following manner. For p = 0 and p = 1,
#' the minimum and maximum of x is used. For other probabilities the quantiles are obtained
#' via \code{q[i] = uniroot(x, cdf - p[i])$root}. Y values are interpolated via
#' \code{\link{approxfun}}.
#'
#' @return Quantiles of cumulative distribution function(s). If the input was a matrix
#' of cumulative distributions functions, a matrix of quantiles is returned.
#'
#' @examples
#' x = seq(-9, 9, .1) # x-grid
#' d = dnorm(x) # density functions
#' p = seq(0, 1, .2) # probabilities of interest
#' cEst = estCdf(d) # estimate cumulative distribution functions
#' qEst = estQdf(p = p, x = x, cdf = cEst) # estimate quantiles
#' plot(x, cEst, bty = 'n', las = 1, type = 'l', ylab = 'Probability') # plot cdf
#' abline(h = p, v = qEst, col = 1:6, lty = 2) # add lines for p and for obtained quantiles
#' points(x = qEst, y = p, pch = 18, col = 1:6, cex = 1.75) # add points for intersections
#'
#' @export
estQdf <- function(p, x, cdf) {
  if (is.matrix(cdf)) {
    out <- apply(cdf, 2, estQdf, p = p, x = x)
    colnames(out) <- colnames(cdf)
    return(out)
  } else if (is.vector(cdf, mode = "numeric")) {
    stopifnot(c(0, 1) %in% p)
    q <- numeric(length = length(p))
    q[1] <- min(x)
    q[length(q)] <- max(x)
    for (i in 2:(length(p) - 1)) {
      f <- stats::approxfun(x, cdf - p[i])
      q[i] <- stats::uniroot(f, range(x))$root
    }
    return(q)
  } else {
    stop(sprintf("Argument cdf must be either a vector or a matrix but it is neither. class(%s) = %s.", 
      deparse(substitute(cdf)), class(cdf)), call. = FALSE)
  }
}


