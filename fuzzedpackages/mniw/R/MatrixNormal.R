#--- matrix normal distribution --------------------------------------------

#' The Matrix-Normal distribution.
#'
#' Density and random sampling for the Matrix-Normal distribution.
#'
#' @name MatrixNormal
#' @aliases dMNorm rMNorm
#' @template param-n
#' @template param-Xpq
#' @template param-Lambda
#' @template param-SigmaR
#' @template param-SigmaC
#' @template param-log
#'
#' @template details-matnorm
#' @example examples/MatrixNormal.R
#'
#' @template return-rdpq

#--- lower-level functions -------------------------------------------------

#' @rdname MatrixNormal
#' @export
dMNorm <- function(X, Lambda, SigmaR, SigmaC, log = FALSE) {
  # get dimensions
  PQ <- .getPQ(X = X, Lambda = Lambda, Sigma = SigmaR, Psi = SigmaC)
  p <- PQ[1]
  q <- PQ[2]
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  X <- .setDims(X, p = p, q = q)
  if(anyNA(X)) stop("Something went wrong.  Please report bug.")
  Lambda <- .setDims(Lambda, p = p, q = q)
  if(anyNA(Lambda)) stop("Lambda and X have incompatible dimensions.")
  SigmaR <- .setDims(SigmaR, p = p)
  if(anyNA(SigmaR)) stop("SigmaR and X have incompatible dimensions.")
  SigmaC <- .setDims(SigmaC, q = q)
  if(anyNA(SigmaC)) stop("SigmaC and X have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, X = X, Lambda = Lambda, Sigma = SigmaR,
             Psi = SigmaC)
  if(length(N) > 2) stop("Arguments have different lengths.")
  ans <- LogDensityMatrixNormal(X, Lambda, SigmaR, SigmaC)
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname MatrixNormal
#' @export
rMNorm <- function(n, Lambda, SigmaR, SigmaC) {
  # get dimensions
  PQ <- .getPQ(Lambda = Lambda, Sigma = SigmaR, Psi = SigmaC)
  p <- PQ[1]
  q <- PQ[2]
  if(any(anyNA(PQ))) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  Lambda <- .setDims(Lambda, p = p, q = q)
  if(anyNA(Lambda)) stop("Something went wrong.  Please report bug.")
  SigmaR <- .setDims(SigmaR, p = p)
  if(anyNA(SigmaR)) stop("SigmaR and Lambda have incompatible dimensions.")
  SigmaC <- .setDims(SigmaC, q = q)
  if(anyNA(SigmaC)) stop("SigmaC and Lambda have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, Lambda = Lambda, Sigma = SigmaR,
             Psi = SigmaC)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n))
    stop("Arguments don't all have length n.")
  X <- GenerateMatrixNormal(n, Lambda, SigmaR, SigmaC)
  if(n > 1) X <- array(X, dim = c(p,q,n))
  X
}
