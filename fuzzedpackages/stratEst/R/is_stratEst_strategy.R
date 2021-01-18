#' Class stratEst.strategy
#' @param x object to be tested.
#' @description Checks if an object is of class \code{stratEst.strategy}.
#' @details Objects of class \code{stratEst.strategy} is returned by the function \code{stratEst.strategy()} of package \code{stratEst}.
#' @export
is.stratEst.strategy <- function( x ){
  inherits(x, "stratEst.strategy")
}
