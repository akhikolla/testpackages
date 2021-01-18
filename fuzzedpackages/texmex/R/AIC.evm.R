#' Information Criteria
#'
#' Compute AIC and (approximate) DIC for \code{evmOpt} objects
#'
#' @aliases AIC.evmSim
#' @param object fit model object
#' @param penalized whether to use the penalized log-likelihood
#' @param nsamp Number of approximate Gaussian sample to use in computing DIC.
#'   Defaults to \code{nsamp=1e3}. Only used when the object has class 'evmOpt'.
#' @param ... other arguments currently ignored
#' @param k numeric, the penalty per parameter to be used; the
#'     default \code{k = 2} is the classical AIC.
#' @param DIC Logical. Whether to compute DIC. Defaults to \code{DIC = TRUE}.
#'   Only applicable to objects of class 'evmSim'.
#' @param WAIC Logical. Whether to compute WAIC. Defaults to \code{WAIC = TRUE}.
#'   Only applicable to objects of class 'evmSim'.
#' @details If the object has class 'evmOpt', \code{nsamp} random draws are
#'   made from the Gaussian distribution with mean and covariance inferred from
#'   the model object. The result will be an approximate DIC. Note that AIC should
#'   not be trusted if priors are not flat. For example, if you use a regularizing
#'   prior on xi, say xi ~ N(0, 0.25), AIC can be misleading and DIC should be
#'   preferred. If the object has class 'evmSim', the actual posterior draws are
#'   used in the computation. Also note that sometimes the optimizer returns
#'   an approximatae covariance that is not postive-semidefinite, in which case
#'   the DIC will be reported as NA.
#' @return The AIC and DIC
#' @seealso \code{\link[stats]{AIC}}
#' @importFrom stats AIC logLik
#' @export
AIC.evmOpt <- function(object, penalized=FALSE, nsamp=1e3, DIC, WAIC, ..., k=2){
  ll <- unclass(logLik(object, penalized=penalized))
  AIC <- -2*ll + k * attr(ll, 'df')
  attr(AIC, "df") <- NULL
  AIC
}

DIC.evmSim <- function(object){
  ll <- object$map$family$log.lik(object$map$data, object$map$th)
  samp <- object$param

  dev <- -2 * apply(samp, 1, ll)
  #dev <- dev[dev < Inf]

  mean(dev) + (mean(dev) + 2 * ll(colMeans(samp)))
}

WAIC.evmSim <- function(object){
  # This follows Richard McElreath's Statistical Rethinking, Section 6.4.2
  # Get density function
  dfun <- object$map$family$density
  samp <- object$param

  # Need matrix with one column for each parameter
  npar <- length(coef(object))

  D <- texmexMakeNewdataD(object$map, newdata=NULL)

  n <- length(object$map$data$y)

  # Get loglik for every row in samp, every value of (phi, xi)
  param <- texmexGetParam(data = D, co = samp)

  lps <- lapply(1:nrow(D[[1]]), # For each observation get matrix of parameters
                function(i, x, p){
                  wh <- lapply(1:length(D),
                               function(j, x, p, i){
                                 rowSums(t(t(p[[j]]) * c(x[[j]][i, ])))
                               }, x=x, p=p, i=i)
                  wh <- do.call("cbind", wh)
                  colnames(wh) <- names(x)
                  wh
                }, x=D, p=param)

  # lps is a list of matrices, one for each observation in the original dataset,
  # each matrix having columns for the parameters (phi and xi for GPD), and
  # with number of rows equal to the number of retained samples from the
  # Markov chains

  ll <- sapply(1:length(lps), function(X){
    dfun(object$map$data$y[X], lps[[X]], object$map, log.d=TRUE)
  })


  lppd <- apply(ll, 2, function(X){
    xm <- max(X)
    s <- sum(exp(X - xm))
    xm + log(s)
  }) - log(n)

  # Get effective number of parameters
  ep <- apply(ll, 2, var)

  -2 * (sum(lppd) - sum(ep))
}


#' @export
AIC.evmSim <- function(object, DIC = TRUE, WAIC = TRUE, ..., k=2){
  aic <- AIC.evmOpt(object$map, ..., k=2)
  aic <- c(AIC = aic)
  if (DIC){
    dic <- DIC.evmSim(object)
    aic <- c(aic, DIC = dic)
  }
  if (WAIC){
    waic <- WAIC.evmSim(object)
    aic <- c(aic, WAIC = waic)
  }

  aic
}

#' Log-likelihood for evmOpt objects
#'
#' Return the log-likelihood or penalized log-likelihood for
#' \code{evmOpt} objects.
#'
#' @param object fit model object
#' @param penalized whether to return the penalized log-likelihood
#' @param ... some methods need more arguments
#' @return an object of class \code{logLik}
#' @seealso \code{\link[stats]{logLik}}
#' @importFrom stats logLik
#' @export
logLik.evmOpt <- function(object, penalized=FALSE, ...){
  ll <- if (penalized) {object$ploglik} else {object$loglik}
  attr(ll, 'df') <- length(coef(object))
  class(ll) <- 'logLik'
  ll
}
