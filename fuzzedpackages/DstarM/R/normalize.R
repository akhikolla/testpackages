#' Normalize two pdfs
#' @param x Probability density function(s) evaluated at grid \code{x}.
#' Input should be either a vector or matrix. If input is a matrix, each column represents a single pdf.
#' @param tt a numeric grid defined in \code{x}.
#' @param props the value each density should integrate to.
#' @examples
#' tt <- seq(0, 9, length.out = 1e4)
#' # 2 poper densities
#' x1 <- cbind(dexp(tt, .5), dexp(tt, 2))
#' # still 2 poper densities
#' x2 <- normalize(10*x1, tt)
#' # 2 densities that integrate to .5
#' x3 <- normalize(x1, tt, props = c(.5, .5))
#' # plot the results
#' matplot(tt, cbind(x1, x2, x3), type = "l", ylab = "density",
#'         col = rep(1:3, each = 2), lty = rep(1:2, 3), las = 1, bty = "n")
#' legend("topright", legend = rep(paste0("x", 1:3), each = 2),
#'        col = rep(1:3, each = 2), lty = rep(1:2, 3), bty = "n")

#' @export
normalize <- function(x, tt, props = NULL) {

  x <- as.matrix(x)
  nc <- dim(x)[2L]
  if (nc != 1L && nc != 2L)
    stop("x must be a matrix with 2 columns.")
  if (any(is.infinite(x), is.na(x), x < 0))
    stop("x contains missing, infinite, or negative values.")
  if (!is.null(props) && any(is.infinite(props), is.na(props), props < 0))
    stop("props contains missing, infinite, or negative values.")

  norm <- diag(nc) / apply(x, 2L, simpson, x = tt)
  if (!is.null(props))
    norm <- norm * props

  return(x %*% norm)
}
