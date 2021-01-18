
#' Random generation from product kernel density
#'
#' @param y         numeric matrix or data.frame.
#' @param n         number of observations. If \code{length(n) > 1},
#'                  the length is taken to be the number required.
#' @param bw        numeric vector of length equal to \code{ncol(y)};
#'                  the smoothing bandwidth to be used. The kernels are
#'                  scaled such that this is the standard deviation of the
#'                  smoothing kernel (see
#'                  \code{\link[stats]{density}} for details). If provided as
#'                  a single value, the same bandwidth is used for each variable.
#' @param weights   numeric vector of length equal to \code{nrow(y)}; must be non-negative.
#' @param adjust    scalar; the bandwidth used is actually \code{adjust*bw}.
#'                  This makes it easy to specify values like 'half the default'
#'                  bandwidth.
#' @param kernel    a character string giving the smoothing kernel to be used.
#'                  This must partially match one of "gaussian", "rectangular",
#'                  "triangular", "epanechnikov", "biweight", "cosine"
#'                  or "optcosine", with default "gaussian", and may be abbreviated.
#' @param shrinked  if \code{TRUE} random generation algorithm preserves mean and
#'                  variances of the individual variables (see \code{\link{ruvk}}).
#'                  Shrinking is applied to each of the variables individually.
#'
#'
#' @details
#'
#' Product kernel density is defined in terms of independent univariate kernels
#'
#' \deqn{
#' \hat{f_H}(\mathbf{x}) = \sum_{i=1}^n w_i \prod_{j=1}^m
#' K_{h_j}( x_i - y_{ij})
#' }{
#' f(x) = sum[i](w[i] * prod[j]( Khj(x[i]-y[i,j]) ))
#' }
#'
#' where \eqn{w} is a vector of weights such that all \eqn{w_i \ge 0}{w[i] \ge 0}
#' and \eqn{\sum_i w_i = 1}{sum(w) = 1} (by default uniform \eqn{1/n} weights are used),
#' \eqn{K_{h_j}}{Khj} is univariate kernel \eqn{K} parametrized by bandwidth \eqn{h_j}{h[j]},
#' where \eqn{\boldsymbol{y}}{y} is a matrix of data points used for estimating the
#' kernel density.
#'
#' For functions estimating kernel densities please check \pkg{KernSmooth},
#' \pkg{ks}, or other packages reviewed by Deng and Wickham (2011).
#'
#' For random generation the algorithm described in \code{\link{kernelboot}} is used.
#' When using \code{shrinked = TRUE}, random noise is drawn from independent, shrinked
#' univariate kernels.
#'
#'
#' @references
#' Deng, H. and Wickham, H. (2011). Density estimation in R.
#' \url{http://vita.had.co.nz/papers/density-estimation.pdf}
#'
#'
#' @examples
#'
#' dat <- mtcars[, c("mpg", "disp")]
#'
#' partmp <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
#'
#' plot(rmvk(5000, dat, shrinked = FALSE), col = "#458B004D", pch = 16,
#'      xlim = c(0, 45), ylim = c(-200, 800),
#'      main = "Product kernel", axes = FALSE)
#' points(dat, pch = 2, lwd = 2, col = "red")
#' axis(1); axis(2)
#'
#' plot(rmvk(5000, dat, shrinked = TRUE), col = "#458B004D", pch = 16,
#'      xlim = c(0, 45), ylim = c(-200, 800),
#'      main = "Product kernel (shrinked)", axes = FALSE)
#' points(dat, pch = 2, lwd = 2, col = "red")
#' axis(1); axis(2)
#'
#' par(partmp)
#'
#' cov(dat)
#' cov(rmvk(5000, dat, shrinked = FALSE))
#' cov(rmvk(5000, dat, shrinked = TRUE))
#'
#' @seealso \code{\link{kernelboot}}
#'
#' @export

rmvk <- function(n, y, bw = sqrt(diag(bw.silv(y))),
                  kernel = c("gaussian", "epanechnikov", "rectangular",
                             "triangular", "biweight", "cosine", "optcosine"),
                  weights = NULL, adjust = 1, shrinked = FALSE) {

  kernel <- match.arg(kernel)
  if (length(n) > 1L) n <- length(n)

  if (is.simple.vector(y)) {
    y <- matrix(y, nrow = 1L)
  } else {
    y <- as.matrix(y)
  }

  if (is.null(weights)) weights <- 1

  if (is.matrix(bw) || is.data.frame(bw)) {
    if (!is.square(bw))
      stop("bw is not a square matrix")
    bw <- diag(bw)
  }
  bw <- bw * adjust[1L]
  if (length(bw) == 1L)
    bw <- rep(bw, ncol(y))

  out <- cpp_rmvk(n, y, bw, weights, kernel, shrinked)
  attr(out, "boot_index") <- NULL
  colnames(out) <- colnames(y)
  out

}
