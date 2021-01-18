#' Calculate upper end point for a fitted extreme value model
#' 
#' Calculate upper end point for fitted extreme value model
#' 
#' 
#' @aliases endPoint endPoint.evmOpt endPoint.evmSim endPoint.evmBoot
#' @usage endPoint(y, verbose=TRUE, .unique=TRUE, ...)
#' 
#' \method{endPoint}{evmOpt}(y,verbose=TRUE, .unique=TRUE, ...)
#' \method{endPoint}{evmSim}(y,verbose=TRUE, .unique=TRUE, ...)
#' @param y Object of class \code{evmOpt} or \code{evmSim}, as returned by
#' \code{\link{evm}}.
#' @param verbose Whether to print output.
#' @param .unique Whether or not to use only unique values of \code{y}.
#' @param ... further arguments to be passed to the \code{\link{signif}}
#' function.
#' @return In cases where the fitted shape parameter is negative, the fitted
#' finite upper endpoint of the extreme value model.
#' @author Janet E. Heffernan
#' @export endPoint
endPoint <- function(y,verbose=TRUE,.unique=TRUE,...){
  UseMethod("endPoint", y)
}

#' @export
endPoint.evmOpt <- function(y, verbose=TRUE,.unique=TRUE,...){

  if(.unique) Unique <- unique else Unique <- identity

  p <- texmexMakeParams(coef(y), y$data$D)
  endpoint <- y$family$endpoint

  negShape <- p[, ncol(p)] < 0

  if(any(negShape)){
    UpperEndPoint <- endpoint(p, y)
    UpperEndPoint[!negShape] <- Inf
    if(verbose){
      o <- Unique(cbind(y$data$D[['xi']], p))
      print(signif(o,...))
    } else {
      invisible(Unique(UpperEndPoint))
    }
  } else {
      Unique(rep(Inf,length(negShape)))
  }
}

#' @export
endPoint.evmSim <- function(y,verbose=TRUE,.unique=TRUE,...){
  endPoint(y$map,verbose=verbose,.unique=.unique,...)
}

#' @export
endPoint.evmBoot <- endPoint.evmSim
