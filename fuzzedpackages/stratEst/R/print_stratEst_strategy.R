#' Print Method for stratEst.strategy
#' @param x An object of class \code{stratEst.strategy}.
#' @param ... Further arguments passed to or from other methods.
#' @export

print.stratEst.strategy <- function( x , ... ){
  x <- round.stratEst.strategy(x,digits=3)
  print.data.frame(x , ... )  #unlist("stratEst.strategy", NextMethod())

}
