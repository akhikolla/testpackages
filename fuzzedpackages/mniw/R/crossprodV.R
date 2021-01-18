#' Matrix cross-product.
#'
#' Vectorized matrix cross-products \code{t(X) V Y} or \code{t(X) V^{-1} Y}.
#'
#' @param X A matrix of size \code{p x q}, or an array of size \code{p x q x n}.
#' @param Y A matrix of size \code{p x r}, or an array of size \code{p x r x n}.  If missing defaults to \code{Y = X}.
#' @param V A matrix of size \code{p x p}, or an array of size \code{p x p x n}.
#' @param inverse Logical; whether or not the inner product should be calculated with \code{V} or \code{V^{-1}}.
#' @return An array of size \code{q x r x n}.
#'
#' @example examples/crossprodV.R
#' @export
crossprodV <- function(X, Y = NULL, V, inverse = FALSE) {
  # dimensions of X and V
  if(is.vector(X)) X <- as.matrix(X)
  p <- dim(X)[1]
  q <- dim(X)[2]
  Xnames <- dimnames(X)[2]
  X <- matrix(X,p)
  if(!all(dim(V)[1:2] == p)) stop("X and V have incompatible dimensions.")
  V <- matrix(V,p)
  if(is.null(Y)) {
    # no Y
    r <- q
    Ynames <- Xnames
    # check lengths
    n <- unique(sort(c(ncol(X)/q, ncol(V)/p)))
    n <- c(1, n[n>1])
    if(length(n) > 2) stop("X and V have different lengths.")
    n <- max(n)
    # C code
    W <- CrossProdVXX(X, V, p, q, inverse)
  } else {
    # dimensions of Y
    if(is.vector(Y)) Y <- as.matrix(Y)
    if(dim(Y)[1] != p) stop("Y and V have incompatible dimensions.")
    Ynames <- dimnames(Y)[2]
    r <- dim(Y)[2]
    Y <- matrix(Y,p)
    # check lengths
    n <- unique(sort(c(ncol(X)/q, ncol(V)/p, ncol(Y)/r)))
    n <- c(1, n[n>1])
    if(length(n) > 2) stop("X, Y, and V have different lengths.")
    n <- max(n)
    # C code
    W <- CrossProdVXY(X, Y, V, p, q, r, inverse)
  }
  W <- array(W, dim = c(q,r,n))
  if(!is.null(Xnames) || !is.null(Ynames)) {
    dimnames(W) <- c(Xnames, Ynames, NULL)
  }
  W
}
