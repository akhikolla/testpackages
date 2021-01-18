# Following simulate.lm, simulate new set of responses from the whole data



#' Simulate from a fitted evm object
#' 
#' Simulate random numbers from a fitted evm object
#' 
#' For \code{simulate.evmSim} and \code{simulate.evmBoot}, the parameters from
#' the Markov chains or bootstrap replicates are randomly permuted prior to
#' each set of simulated responses being computed. In this way, reusing the
#' same set of values is avoided.
#' 
#' @param object A fitted evm object having class 'evmOpt', 'evmSim' or
#' 'evmBoot'.
#' @param nsim The number of simulations to perform. Defaults to \code{nsim=1}.
#' A single simulation involves simulating a new set of responses from the data
#' that was provided to \code{evm} (after thresholding if thresholding is
#' performed.)
#' @param seed An integer to be passed to \code{set.seed}. Defaults to
#' \code{seed=NULL}.
#' @param param Parameters to use in the random number generator. Defaults to
#' \code{param=NULL} in which case the parameters from the fitted model are
#' used.  For \code{simulate.evmSim} and \code{simulate.evmBoot}, this argument
#' is not available and the simulated parameters or replicates are used.
#' @param ... Unused.
#' @return If \code{nsim=1}, a vector or random numbers simulated from the
#' fitted model object.  If \code{nsim > 1}, a matrix with each column being a
#' set of simulated responses.
#' @author Paul Metcalfe, Harry Southworth
#' @seealso \code{\link{evm}}
#' @keywords models
#' @examples
#' 
#' mod <- evm(rain, qu=.95)
#' hist(simulate(mod, 100))
#' 
#' @export
simulate.evmOpt <- function(object, nsim=1, seed=NULL, param=NULL, ...){
  if (!is.null(seed)) set.seed(seed)
  if (is.null(param)){
    param <- predict(object, type="lp", unique.=FALSE)$obj
    param <- param$link[, colnames(param$link) %in% names(param$family$param)]
    if (nsim > 1)
      param <- do.call("rbind", rep(list(param), nsim))
  }

  res <- unname(object$family$rng(nrow(param), param, object))
  if (nsim > 1)
    matrix(res, byrow=FALSE, ncol=nsim)
  else
    res
}

#' @rdname simulate.evmOpt
#' @export
simulate.evmSim <- function(object, nsim=1, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)

  res <- list()
  for (i in 1:nsim){ # Do this here to stop reuse of same parameters
    # permute parameters to avoid reusing the same ones when nsim > 1
    object$param <- object$param[sample(nrow(object$param)), ]
    
    param <- predict(object, type="lp", unique.=FALSE)$obj
    param <- param$link[, colnames(param$link) %in% names(param$family$param)]
    
    res[[i]] <- simulate(object$map, nsim=1, param=param, seed=NULL)
  }
  do.call("cbind", res)
}

#' @rdname simulate.evmOpt
#' @export
simulate.evmBoot <- function(object, nsim=1, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)

  res <- list()
  for (i in 1:nsim){
    # permute parameters to avoid reusing the same ones when nsim > 1
    object$replicates <- object$replicates[sample(nrow(object$replicates)), ]
    
    param <- predict(object, type="lp", unique.=FALSE)$obj
    param <- param$link[, colnames(param$link) %in% names(param$family$param)]
    
    res[[i]] <- simulate(object$map, nsim=1, param=param, seed=NULL)
  }
  do.call("cbind", res)
}
