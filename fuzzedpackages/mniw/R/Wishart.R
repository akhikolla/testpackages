#---------------------------------------------------------------------------
#
# Wishart and Inverse-Wishart distributions
#
#---------------------------------------------------------------------------

#' Wishart and Inverse-Wishart distributions.
#'
#' Densities and random sampling for the Wishart and Inverse-Wishart distributions.
#'
#' @name Wishart
#' @aliases dwish rwish diwish riwish dwishart rwishart
#' @template param-Xqq
#' @template param-n
#' @template param-Psi
#' @template param-nu
#' @param inverse Logical; whether or not to use the Inverse-Wishart distribution.
#' @template param-log
#' @template details-wishart
#' @details \code{dwish} and \code{diwish} are convenience wrappers for \code{dwishart}, and similarly \code{rwish} and \code{riwish} are wrappers for \code{rwishart}.
#'
#' @example examples/Wishart.R
#' @template return-rdqq

#--- convenience wrappers --------------------------------------------------

# \deqn{
# f(\boldsymbol{X} \mid \boldsymbol{\Psi}, \nu) = \frac{|\boldsymbol{X}|^{(\nu-q+1)/2}\exp\big\{-\frac 1 2 \textrm{trace}(\boldsymbol{\Psi}^{-1}\boldsymbol{X})\big\}}{2^{\nu q/2}|\boldsymbol{\Psi}|^{\nu/2} \Gamma_q(\frac \nu 2)}
# }


# wishart density
#' @rdname Wishart
#' @export
dwish <- function(X, Psi, nu, log = FALSE) {
  dwishart(X, Psi, nu, inverse = FALSE, log)
}

# wishart simulation
#' @rdname Wishart
#' @export
rwish <- function(n, Psi, nu) {
  rwishart(n, Psi, nu, inverse = FALSE)
}

# inverse wishart density
#' @rdname Wishart
#' @export
diwish <- function(X, Psi, nu, log = FALSE) {
  dwishart(X, Psi, nu, inverse = TRUE, log)
}

# inverse wishart simulation
#' @rdname Wishart
#' @export
riwish <- function(n, Psi, nu) {
  rwishart(n, Psi, nu, inverse = TRUE)
}

#--- lower level functions -------------------------------------------------

# density of wishart and inverse wishart

#' @rdname Wishart
#' @export
dwishart <- function(X, Psi, nu, inverse = FALSE, log = FALSE) {
  # get dimensions
  PQ <- .getPQ(X = X, Psi = Psi)
  if(is.na(PQ[2])) stop("Undetermined problem dimensions.")
  q <- PQ[2]
  # format arguments
  X <- .setDims(X, q = q)
  Psi <- .setDims(Psi, q = q)
  nu <- c(nu)
  if(anyNA(X)) stop("Something went wrong.  Please report bug.")
  if(anyNA(Psi)) stop("Psi and X have incompatible dimensions.")
  if(length(X) == 1) X <- as.matrix(X)
  if(length(Psi) == 1) Psi <- as.matrix(Psi)
  # check lengths
  N <- .getN(q = q, X = X, Psi = Psi, nu = nu)
  ## check X
  #q <- dim(X)[1]
  #if(dim(X)[2] != q) stop("X has non-square dimensions.")
  #X <- matrix(X,q)
  ## check Psi
  #if(missing(Psi)) Psi <- diag(q)
  #if(!all(dim(Psi)[1:2] == q)) stop("X and Psi have incompatible dimensions.")
  #Psi <- matrix(Psi,q)
  ## check nu
  #nu <- c(nu)
  ## check lengths
  #N <- unique(sort(c(ncol(X)/q, ncol(Psi)/q, length(nu))))
  #N <- c(1, N[N>1])
  if(length(N) > 2) stop("Arguments have different lengths.")
  ans <- LogDensityWishart(X, Psi, nu, inverse)
  if(!log) ans <- exp(ans)
  ans
}

# random sampling for wishart and inverse wishart
#' @rdname Wishart
#' @export
rwishart <- function(n, Psi, nu, inverse = FALSE) {
  # get problem dimensions
  PQ <- .getPQ(Psi = Psi)
  if(is.na(PQ[2])) stop("Undetermined problem dimensions.")
  q <- PQ[2]
  # format arguments
  Psi <- .setDims(Psi, q = q)
  nu <- c(nu)
  # check lengths
  N <- .getN(q = q, Psi = Psi, nu = nu)
  ## determine the dimensions
  #q <- dim(Psi)[1]
  #if(dim(Psi)[2] != q) stop("Psi has non-square dimensions.")
  #Psi <- matrix(Psi,q)
  #nu <- c(nu)
  ## determine lengths
  #N <- sort(unique(c(ncol(Psi)/q, length(nu))))
  #N <- c(1,N[N>1])
  if(length(N) > 2 || (length(N) == 2 && N[2] != n))
    stop("Arguments don't all have length n.")
  X <- GenerateWishart(n, Psi, nu, inverse)
  if(n > 1) X <- array(X, dim = c(q,q,n))
  X
}

