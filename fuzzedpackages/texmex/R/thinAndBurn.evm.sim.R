#' Process Metropolis output from extreme value model fitting to discard
#' unwanted observations.
#'
#' Process observations from Metropolis fitting of extreme value models, to
#' thin the output and discard observations from burn-in period.
#'
#'
#' @aliases thinAndBurn thinAndBurn.evmSim
#' @usage \method{thinAndBurn}{evmSim}(object, burn, thin)
#' @param object Object of class 'evmSim' as returned by \code{evm} called with
#' \code{method="simulate"}.
#' @param thin \code{thin} or its reciprocal must be a positive integer.  If
#' integer valued, this specifies the frequency of observations from the
#' simulated Markov Chain which will be retained.  If specified as a
#' proportion, this is the proportion of values which will be retained. For no
#' thinning use \code{thin=1}.
#' @param burn The number of observations from the simulated Markov Chain to be
#' discarded as burn-in. Must be a non-negative integer, for no burn-in use
#' \code{burn=0}.
#' @return Object of class \code{evmSim}.  See Value returned by
#' \code{\link{evm}} using \code{method = "simulate"} for details.
#'
#' Note that the original chain is not discarded when this function is called:
#' \code{thinAndBurn} can be called recursively on the original object with
#' different values of \code{burn} and \code{thin} without the object getting
#' progressively smaller!
#' @author Harry Southworth, Janet E. Heffernan
#' @seealso \code{\link{evm}}
#' @examples
#'
#'   x <- rnorm(1000)
#'   # For the values of burn and thin below, we should do many more iterations.
#'   # The number of iterations is kept low here due to the run time allowed
#'   # by CRAN.
#'   mod <- evm(x, qu=.7, method="sim", iter=11000)
#'   mod
#'   par(mfrow=c(3, 2))
#'   plot(mod)
#'   mod1 <- thinAndBurn(mod,burn=1000, thin=5)
#'   plot(mod1)
#'
#' @export thinAndBurn
thinAndBurn <- function (object, burn, thin){
  UseMethod("thinAndBurn")
}

#' @export
thinAndBurn.evmSim <- function(object, burn, thin){

  if(missing(burn)){
    burn <- object$burn
  } else {
    object$burn <- burn
  }
  if(missing(thin)){
    thin <- object$thin
  } else {
    object$thin <- thin
  }
  if(is.null(object$thin)){
    stop("thin or its reciprocal must be a positive integer, for no thinning use thin=1")
  }
  if(is.null(object$burn)){
    stop("burn must be a non-negative integer, for no burn in use burn=0")
  }

  if (thin < 1) thin <- 1 / thin
  if (thin %% 1 > 10^(-6)) stop("thin, or its reciprocal, should be an integer")

  if (burn > dim(object$chains[[1]])[1]) stop("burn-in is longer than the whole chain")

  if (burn > 0){
     object$param <- lapply(object$chains, function(X){
       X[-(1:burn), ] # Remove burn-in
     })
  } else {
     object$param <- object$chains
  }

  wh <- 1:nrow(object$param[[1]]) %% thin == 0
  object$param <- lapply(object$param, function(X) X[wh ,])

  object$param <- do.call("rbind", object$param)

  invisible(object)
}

