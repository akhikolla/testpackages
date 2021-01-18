#' Bayesian inference for a random-effects regression model.
#'
#' Gibbs sampler for posterior distribution of parameters and hyperparameters of a multivariate normal random-effects linear regression model called RxNormLM (see \strong{Details}).
#'
#' @param nsamples number of posterior samples to draw.
#' @param Y \code{N x q} matrix of responses.
#' @param V Either a \code{q x q} variance matrix or an \code{q x q x N} array of such matrices.
#' @param X \code{N x p} matrix of covariates.
#' @param prior parameters of the prior MNIW distribution on the hyperparameters (see \strong{Details}).
#' @param burn integer number of burn-in samples, or fraction of \code{nsamples} to prepend as burn-in.
#' @param init (optional) list with elements \code{Beta}, \code{Sigma}, and \code{Mu} providing the initial values for these.  Default values are \code{Beta = matrix(0, p, q)}, \code{Sigma = diag(q)}, and \code{Mu = Y}.
#' @param updateHyp,storeHyp logical. Whether or not to update/store the hyperparameter draws.
#' @param updateRX,storeRX logical. Whether or not to update/store the random-effects draws.
#' @details The RxNormLM model is given by
#' \deqn{
#' y_i \mid \mu_i \sim_iid N(\mu_i, V_i)
#' }
#' \deqn{
#' \mu_i \mid \beta, \Sigma ~sim_ind N(x_i' \beta, \Sigma)
#' }
#' \deqn{
#' \beta, \Sigma ~ MNIW(\Lambda, \Omega^{-1}, \Psi, \nu),
#' }
#' where \eqn{y_i} and \eqn{\mu_i} are response and random-effects vectors of length \eqn{q}, \eqn{x_i} are covariate vectors of length \eqn{p}, and \eqn{(\beta, \Sigma)} are hyperparameter matrices of size \eqn{p \times q} and \eqn{q \times q}.
#'
#' The MNIW prior distribution is given by a list with elements \code{Lambda}, \code{Omega}, \code{Psi}, and \code{nu}.  If any of these is \code{NULL} or missing, the default value is 0.  Note that \code{Omega == 0} gives a Lebesgue prior to \eqn{\beta}.
#' @return A list with (potential) elements:
#' \describe{
#'   \item{\code{Beta}}{An \code{p x q x nsamples} array of regression coefficient iterations (if \code{storeHyp == TRUE})}
#'   \item{\code{Sigma}}{An \code{q x q x nsamples} array of regression variance matrices (if \code{storeHyp == TRUE})}
#'   \item{\code{Mu}}{An \code{n x q x nsamples} array of random effects (if \code{storeRX == TRUE})}
#' }
#'
#' @example examples/Hierarchical.R
#' @export
RxNormLM <- function(nsamples, Y, V, X, prior = NULL, init, burn,
                      updateHyp = TRUE, storeHyp = TRUE,
                      updateRX = TRUE, storeRX = FALSE) {
  # argument check
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  if(nrow(X) != n) stop("Y and X have incompatible dimensions.")
  if(is.matrix(V)) V <- array(c(V), dim = c(dim(V), n))
  V <- .setDims(V, q, q)
  if(anyNA(V) || !identical(dim(V), c(q,q*n))) {
    stop("Y and V have incompatible dimensions.")
  }
  # prior
  if(is.null(prior)) prior <- list()
  Lambda <- .setDims(prior$Lambda, p, q)
  if(anyNA(Lambda)) stop("Lambda and (Y, X) have incompatible dimensions.")
  Omega <- .setDims(prior$Omega, p, p)
  if(anyNA(Omega)) stop("Omega and X have incompatible dimensions.")
  Psi <- .setDims(prior$Psi, q, q)
  if(anyNA(Psi)) stop("Psi and Y have incompatible dimensions.")
  nu <- if(is.null(prior$nu)) 0 else prior$nu
  if(length(nu) != 1) stop("nu must be a scalar.")
  # initial values
  if(is.null(init)) init <- list()
  Beta0 <- .setDims(init$Beta, p, q)
  if(anyNA(Beta0)) {
    stop("init$Beta and (Y, X) have incompatible dimensions.")
  }
  Sigma0 <- .setDims(init$Sigma, q)
  if(anyNA(Sigma0)) {
    stop("init$Sigma and Y have incompatible dimensions.")
  }
  Mu0 <- if(is.null(init$Mu)) Y else init$Mu
  if(!identical(dim(Mu0), c(n, q))) {
    stop("init$Mu and Y have incompatible dimensions.")
  }
  # burn-in
  if (missing(burn)) burn <- max(0.1, 1000)
  if (burn < 1) burn <- nsamples * burn
  burn <- floor(burn)
  # don't store if don't update
  if(!updateHyp) storeHyp <- FALSE
  if(!updateRX) storeRX <- FALSE
  # MCMC
  post <- HierUneqVModelGibbs(nSamples = as.integer(nsamples),
                              nBurn = as.integer(burn),
                              Y = Y, X = X, V = V,
                              Lambda = Lambda, Omega = Omega,
                              Psi = Psi, nu = nu,
                              Beta0 = Beta0, iSigma0 = .solveV(Sigma0),
                              Mu0 = Mu0,
                              updateBetaSigma = as.logical(updateHyp),
                              updateMu = as.logical(updateRX),
                              storeBetaSigma = as.logical(storeHyp),
                              storeMu = as.logical(storeRX))
  # output format
  out <- NULL
  if(storeHyp) {
    out <- list(Beta = array(post$Beta, dim = c(p, q, nsamples)),
                Sigma = array(post$Sigma, dim = c(q, q, nsamples)))
  }
  if(storeRX) {
    out <- c(out,
             list(Mu = aperm(array(post$Mu, dim = c(q, n, nsamples)),
                             perm = c(2,1,3))))
  }
  out
}

# solve for variance matrices
.solveV <- function(V, x) {
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
}

