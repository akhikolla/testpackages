#' Class stratEst.check
#' @param x object to be tested.
#' @description Checks if an object is of class \code{stratEst.check}.
#' @details Objects of class \code{stratEst.check} are returned by the function \code{stratEst.check()} of package \code{stratEst}.
#' @export
is.stratEst.check <- function( x ){
  inherits(x, "stratEst.check")
}
