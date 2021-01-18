#' Class stratEst.data
#' @param x object to be tested.
#' @description Checks if an object is of class \code{stratEst.data}.
#' @details Objects of class \code{stratEst.data} are returned by the functions \code{stratEst.data()} and \code{stratEst.simulate()} of package \code{stratEst}.
#' @export
is.stratEst.data <- function( x ){
  inherits(x, "stratEst.data")
}
