#' Sparse projections using the square-root Lasso
#'
#' Regresses each column of \code{x} against all others in turn, using the
#' square-root Lasso, and outputs residuals from the regressions. Thus it
#' outputs a form of sparse projection of each column on all others.
#' Alternatively, given two matrices \code{x_null} and \code{x_alt}, it
#' regresses each column of \code{x_null} on \code{x_alt} in a similar fashion.
#'
#' @param x Matrix with each row an observation vector. Need not be supplied if
#'   \code{x_alt} and \code{x_null} are given.
#' @param x_null Matrix whose columns are to be regressed on to \code{x_alt}.
#' @param x_alt Matrix which the columns of \code{x_null} are regressed on to.
#'   Must be specified if \code{x_null} is given.
#' @param mc.cores The number of cores to use. Will always be 1 in Windows.
#' @param rescale Should the columns of the output be rescaled to have l_2-norm
#'   the square-root of the number of observations? Default is \code{FALSE}.
#' @param ... Additional arguments to be passed to \code{sqrt_lasso}.
#' @return A matrix where each column gives the residuals.
#' @references A. Belloni, V. Chernozhukov, and L. Wang. (2011)
#'   \emph{Square-root lasso: pivotal recovery of sparse signals via conic
#'   programming. Biometrika, 98(4):791-806.}
#'   \url{http://biomet.oxfordjournals.org/content/98/4/791.refs} T. Sun and
#'   C.-H. Zhang. (2012) \emph{Scaled sparse linear regression. Biometrika,
#'   99(4):879-898.}
#'   \url{http://biomet.oxfordjournals.org/content/early/2012/09/24/biomet.ass043.short}
#'    T. Sun and C.-H. Zhang. (2013) \emph{Sparse matrix inversion with scaled
#'   lasso. The Journal of Machine Learning Research, 14(1):3385-3418.}
#'   \url{www.jmlr.org/papers/volume14/sun13a/sun13a.pdf}
#' @seealso \code{\link{sqrt_lasso}} and \code{\link{RPtest_single}}.
#' @examples
#' x <- matrix(rnorm(50*100), 50, 100)
#' out <- sparse_proj(x)
#' @export
sparse_proj <- function(x, x_null, x_alt, mc.cores=1L, rescale=FALSE, ...) {
  if (missing(x_null)) {
    # Check x
    if (!is.matrix(x)) stop("x should be a matrix with two or more columns")
    np <- dim(x)
    if (is.null(np) | (np[2] < 2L))
      stop("x should be a matrix with at least two columns")

    n <- nrow(x)
    p <- ncol(x)
    x <- scale(x)*sqrt(n/(n-1))
    rep_fun <- function(var) {
      resid_lasso(x[, -var], x[, var], rescale=rescale, ...)
    }
    out <- parallel::mclapply(1:p, rep_fun, mc.cores=mc.cores)
    out <- matrix(unlist(out), nrow=n, ncol=p)
    return(out)
  } else {
    if (missing(x_alt)) stop("x_alt must be specified if x_null is given")
    # Check x_alt and x_null
    if (!(is.matrix(x_null) && is.matrix(x_alt))) stop("x_null and x_alt must be matrices with at least one column")
    n <- nrow(x_alt); p <- ncol(x_alt)
    rep_fun <- function(var) {
      resid_lasso(x_null, x_alt[, var], rescale=rescale, ...)
    }
    out <- parallel::mclapply(1:ncol(x_alt), rep_fun, mc.cores=mc.cores)
    out <- matrix(unlist(out), nrow=n, ncol=p)
    return(out)
  }
}

# mc.cores and resid_type are unused arguments included simply to match with sparse_proj
non_sparse_proj <- function(x_null, x_alt, mc.cores = NULL, resid_type=NULL) {
  u <- svd(x_null, nv=0)$u
  x_alt <- x_alt - u %*% t(u) %*% x_alt
  return(x_alt)
}
