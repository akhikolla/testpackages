#--- MNIW distribution -----------------------------------------------------

#' Generate samples from the Matrix-Normal Inverse-Wishart distribution.
#'
#' @name MNIW
#' @param n number of samples.
#' @param Lambda A mean matrix of size \code{p x q} or an array of size \code{p x q x n}.  Defaults to matrix of zeros when missing.
#' @param Sigma A row-wise variance or precision matrix of size \code{p x p}, or an array of size \code{p x p x n}.  Defaults to the identity matrix when missing.
#' @param Omega A between-row precision matrix of size \code{p x p}, or an array of size \code{p x p x n}.  Defaults to the identity matrix when missing.
#' @param Psi A scale matrix of size \code{q x q}, or an array of size \code{q x q x n}.  Defaults to identity matrix when missing.
#' @param nu Scalar degrees-of-freedom parameter.
#' @param prec Logical; whether or not \code{Sigma} is on the variance or precision scale.
#'
#' @return A list with elements:
#' \describe{
#' \item{\code{X}}{Array of size \code{p x q x n} random samples from the Matrix-Normal component (see \strong{Details}).}
#' \item{\code{V}}{Array of size \code{q x q x n} of random samples from the Inverse-Wishart component.}
#' }
#' @template details-mniw
#' @details \code{rmniw} is a convenience wrapper to \code{rMNIW(Sigma = Omega, prec = TRUE)}, for the common situation in Bayesian inference with conjugate priors when between-row variances are naturally parametrized on the precision scale.
#' @example examples/MNIW.R
#' @export
rMNIW <- function(n, Lambda, Sigma, Psi, nu, prec = FALSE) {
  # get dimensions
  PQ <- .getPQ(Lambda = Lambda, Sigma = Sigma, Psi = Psi)
  p <- PQ[1]
  q <- PQ[2]
  if(anyNA(PQ)) {
    stop("Problem dimensions are undetermined (too many missing inputs).")
  }
  # format arguments
  Lambda <- .setDims(Lambda, p = p, q = q)
  if(anyNA(Lambda)) stop("Something went wrong.  Please report bug.")
  Sigma <- .setDims(Sigma, p = p)
  if(anyNA(Sigma)) stop("Sigma and Lambda have incompatible dimensions.")
  Psi <- .setDims(Psi, q = q)
  if(anyNA(Psi)) stop("Psi and Lambda have incompatible dimensions.")
  nu <- c(nu)
  # check lengths
  N <- .getN(p = p, q = q, Lambda = Lambda, Sigma = Sigma,
             Psi = Psi, nu = nu)
  if(length(N) > 2 || (length(N) == 2 && N[2] != n)) {
    stop("Arguments don't all have length n.")
  }
  XV <- GenerateMatrixNIW(n, Lambda, Sigma, Psi, nu, prec)
  if(n > 1) {
    XV$X <- array(XV$X, dim = c(p,q,n))
    XV$V <- array(XV$V, dim = c(q,q,n))
  }
  XV
}

#' @rdname MNIW
#' @export
rmniw <- function(n, Lambda, Omega, Psi, nu) {
  rMNIW(n = n, Lambda = Lambda, Sigma = Omega,
        Psi = Psi, nu = nu, prec = TRUE)
}
