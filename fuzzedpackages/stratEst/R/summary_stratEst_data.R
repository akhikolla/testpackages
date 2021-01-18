#' Method dispatch for Generic Function summary
#' @param object An object to be summarized.
#' @param ... additional arguments affecting the result.
#' @export

summary.stratEst.data <- function( object , ...){

  c("stratEst.data", NextMethod())

}
