#--- matrix-t distribution ------------------------------------------------

#' The Matrix-t distribution.
#'
#' Density and sampling for the Matrix-t distribution.
#'
#' @name MatrixT
#' @aliases dMT
#' @template param-n
#' @template param-Xpq
#' @template param-Lambda
#' @template param-SigmaR
#' @template param-SigmaC
#' @template param-nu
#' @template param-log
#'
#' @template return-rdpq
#'
#' @template details-matrixt

#--- lower-level functions ------------------------------------------------

#' @rdname MatrixT
#' @export
dMT <- function(X, Lambda, SigmaR, SigmaC, nu, log = FALSE) {
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
  nu <- c(nu)
  # check lengths
  N <- .getN(p = p, q = q, X = X, Lambda = Lambda, Sigma = SigmaR,
             Psi = SigmaC, nu = nu)
  if(length(N) > 2) stop("Arguments have different lengths.")
  ans <- LogDensityMatrixT(X, Lambda, SigmaR, SigmaC, nu)
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname MatrixT
#' @export
rMT <- function(n, Lambda, SigmaR, SigmaC, nu) {
  # get dimensions
  PQ <- .getPQ(Lambda = Lambda, Sigma = SigmaR, Psi = SigmaC)
  p <- PQ[1]
  q <- PQ[2]
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  Lambda <- .setDims(Lambda, p = p, q = q)
  if(anyNA(Lambda)) stop("Something went wrong.  Please report bug.")
  SigmaR <- .setDims(SigmaR, p = p)
  if(anyNA(SigmaR)) stop("SigmaR and Lambda have incompatible dimensions.")
  SigmaC <- .setDims(SigmaC, q = q)
  if(anyNA(SigmaC)) stop("SigmaC and Lambda have incompatible dimensions.")
  nu <- c(nu)
  # check lengths
  N <- .getN(p = p, q = q, Lambda = Lambda, Sigma = SigmaR,
             Psi = SigmaC, nu = nu)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n)) {
    stop("Arguments don't all have length n.")
  }
  X <- GenerateMatrixT(n, Lambda, SigmaR, SigmaC, nu)
  if(n > 1) X <- array(X, dim = c(p,q,n))
  X
}
