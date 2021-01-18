#' Estimate cumulative distribution for D*M models
#'
#' @param x Any density function to calculate a cumulative distribution for.
#' The code is designed for input of class \code{DstarM} but other input is also
#' accepted. Other input can be either a matrix where columns represent densities
#' or a single vector representing a density.
#'
#' @details Cumulative distributions functions are calculated by: \code{cumsum(x) / sum(x)}.
#' This method works well enough for our purposes. The example below shows that the
#' \code{\link{ecdf}} functions seems to work slightly better. However, this estimates a
#' cdf from raw data and does not transform a pdf into a cdf and is therefore not useful
#' for D*M models.
#'
#' @return Cumulative density function(s). If the input was a matrix,
#' a matrix of cumulative density functions is returned.
#'
#' @examples
#' x = rnorm(1000)
#' xx = seq(-5, 5, .1)
#' approx1 = stats::ecdf(x)(xx)
#' approx2 = estCdf(dnorm(xx, mean(x), sd(x)))
#' trueCdf = pnorm(xx)
#' matplot(xx, cbind(trueCdf, approx1, approx2), type = c('l', 'p', 'p'),
#'         lty = 1, col = 1:3, pch = 1, bty = 'n', las = 1, ylab = 'Prob')
#' legend('topleft', legend = c('True Cdf', 'Stats Estatimation', 'DstarM Estimation'),
#'        col = 1:3, lty = c(1, NA, NA), pch = c(NA, 1, 1), bty = 'n')
#'
#' @export
estCdf <- function(x) {
  if (is.DstarM.fitObs(x)) {
    x <- x$obsNorm
  }
  if (is.matrix(x)) {
    out <- apply(x, 2, estCdf)
    colnames(out) <- colnames(x)
    return(out)
  } else {
    return(cumsum(x)/sum(x))
  }
}


