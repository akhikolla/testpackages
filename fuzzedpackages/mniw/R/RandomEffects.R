#--- random effects models -------------------------------------------------

#' Conditional sampling for Multivariate-Normal Random-Effects model.
#'
#' Sample from the conditional parameter distribution given the data and hyperparameters of the Multivariate-Normal Random-Effects (mNormRE) model (see \strong{Details}).
#'
#' @template param-n
#' @param x Data observations.  Either a vector of length \code{q} or a \code{n x q} matrix.  In the latter case each row is a different vector.
#' @param V Observation variances.  Either a matrix of size \code{q x q} or a \code{q x q x n} array.
#' @param lambda Prior means.  Either a vector of length \code{q} or an \code{n x q} matrix.  In the latter case each row is a different mean.  Defaults to zeros.
#' @param Sigma Prior variances.  Either a matrix of size \code{q x q} or a \code{q x q x n} array.  Defaults to identity matrix.
#' @details Consider the hierarchical multivariate normal model
#' \deqn{
#' \arraycolsep=1.4pt
#' \begin{array}{rcl}
#' \boldsymbol{\mu} & \sim & \mathcal{N}(\boldsymbol{\lambda}, \boldsymbol{\Sigma}) \\
#' \boldsymbol{x} \mid \boldsymbol{\mu} & \sim & \mathcal{N}(\boldsymbol{\mu}, \boldsymbol{V}).
#' \end{array}
#' }{
#' \mu ~ N(\lambda, \Sigma)
#' }
#' \deqn{
#' \vspace{-2em}
#' }{
#' x | \mu ~ N(\mu, V).
#' }
#' The Multivariate-Normal Random-Effects model \eqn{\boldsymbol{\mu} \sim \textrm{RxNorm}(\boldsymbol{x}, \boldsymbol{V}, \boldsymbol{\lambda}, \boldsymbol{\Sigma})}{\mu ~ RxNorm(x, V, \lambda, \Sigma)} on the random vector \eqn{\boldsymbol{\mu}_q}{\mu_q} is defined as the posterior distribution \eqn{p(\boldsymbol{\mu} \mid \boldsymbol{x}, \boldsymbol{\lambda}, \boldsymbol{\Sigma})}{p(\mu | x, \lambda, \Sigma)}.  This distribution is multivariate normal; for the mathematical specification of its parameters please see \code{vignette("mniw-distributions", package = "mniw")}.
#'
#' @example examples/RandomEffects.R
#' @export
rRxNorm <- function(n, x, V, lambda, Sigma) {
  # convert to MN format
  x <- .vec2mn(x)
  if(!missing(lambda)) lambda <- .vec2mn(lambda)
  # problem dimensions
  PQ1 <- .getPQ(X = x, Lambda = lambda, Sigma = V)
  PQ2 <- .getPQ(X = x, Sigma = Sigma)
  p <- PQ1[1]
  q <- PQ1[2]
  if(q != 1) stop("x and lambda must be vectors or matrices.")
  # format arguments
  x <- .setDims(x, p = p, q = q)
  lambda <- .setDims(lambda, p = p, q = q)
  if(anyNA(lambda)) stop("lambda and x have incompatible dimensions.")
  V <- .setDims(V, p = p)
  if(anyNA(V)) stop("V and x have incompatible dimensions.")
  Sigma <- .setDims(Sigma, p = p)
  if(anyNA(Sigma)) stop("Sigma and x have incompatible dimensions.")
  # check lengths
  N1 <- .getN(p = p, q = q, X = x, Lambda = lambda, Sigma = V)
  N2 <- .getN(p = p, q = q, X = x, Sigma = Sigma)
  N <- unique(sort(c(N1, N2)))
  N <- c(1, N[N>1])
  if(length(N) > 2 || (length(N) == 2 && N[2] != n)) {
    stop("Arguments don't all have length n.")
  }
  Mu <- GenerateRandomEffectsNormal(n, x, V, lambda, Sigma)
  if(n > 1) {
    Mu <- t(Mu)
  } else {
    Mu <- c(Mu)
  }
  Mu
}
