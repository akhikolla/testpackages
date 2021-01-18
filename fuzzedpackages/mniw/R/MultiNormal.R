#--- multivariate normal distribution --------------------------------------

#' The Multivariate Normal distribution.
#'
#' Density and random sampling for the Multivariate Normal distribution.
#'
#' @name MultiNormal
#' @aliases dmNorm rmNorm
#' @param x Argument to the density function.  A vector of length \code{q} or an \code{n x q} matrix.
#' @template param-n
#' @param mu Mean vector(s).  Either a vector of length \code{q} or an \code{n x q} matrix.  If missing defaults to a vector of zeros.
#' @param Sigma Covariance matrix or matrices.  Either a \code{q x q} matrix or a \code{q x q x n} array.  If missing defaults to the identity matrix.
#' @template param-log
#' @return A vector for densities, or a \code{n x q} matrix for random sampling.
#' @example examples/MultiNormal.R
#' @rdname MultiNormal
#' @export
dmNorm <- function(x, mu, Sigma, log = FALSE) {
  # get dimensions
  # first convert to appropriate MN format
  x <- .vec2mn(x)
  if(!missing(mu)) mu <- .vec2mn(mu)
  PQ <- .getPQ(X = x, Lambda = mu, Sigma = Sigma)
  p <- PQ[1]
  q <- PQ[2]
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  if(q != 1) stop("x and mu must be vectors or matrices.")
  # format arguments
  x <- .setDims(x, p = p, q = q)
  if(anyNA(x)) stop("Something went wrong.  Please report bug.")
  mu <- .setDims(mu, p = p, q = q)
  if(anyNA(mu)) stop("mu and x have incompatible dimensions.")
  Sigma <- .setDims(Sigma, p = p)
  if(anyNA(Sigma)) stop("Sigma and x have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, X = x, Lambda = mu, Sigma = Sigma)
  if(length(N) > 2) stop("Arguments have different lengths.")
  x <- matrix(x, nrow = p) # format for mN
  mu <- matrix(mu, nrow = p) # format for mN
  ans <- LogDensityMultivariateNormal(x, mu, Sigma)
  if(!log) ans <- exp(ans)
  ans
}

#' @rdname MultiNormal
#' @export
rmNorm <- function(n, mu, Sigma) {
  # get dimensions
  # first convert to appropriate MN format
  ## if(missing(mu) && !missing(Sigma)) {
  ##   mu <- rep(0, )
  ## }
  if(!missing(mu)) {
    mu <- .vec2mn(mu)
  }
  PQ <- .getPQ(Lambda = mu, Sigma = Sigma)
  if(is.na(PQ[2])) PQ[2] <- 1 # if mu is missing, still know that q = 1
  p <- PQ[1]
  q <- PQ[2]
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined: must provide mu or Sigma.")
  }
  if(q != 1) stop("mu must be a vector or matrix.")
  # format arguments
  mu <- .setDims(mu, p = p, q = q)
  if(anyNA(mu)) stop("mu and Sigma have incompatible dimensions.")
  Sigma <- .setDims(Sigma, p = p)
  if(anyNA(Sigma)) stop("Sigma and mu have incompatible dimensions.")
  # check lengths
  N <- .getN(p = p, q = q, Lambda = mu, Sigma = Sigma)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n)) {
    stop("Arguments don't all have length n.")
  }
  mu <- matrix(mu, nrow = p) # format for mN
  X <- GenerateMultivariateNormal(n, mu, Sigma)
  if(n > 1) {
    X <- t(X)
  } else {
    X <- c(X)
  }
  X
}
