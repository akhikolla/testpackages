#' Round Method for stratEst.strategy
#' @param x An object of class \code{stratEst.strategy}.
#' @param digits Further arguments passed to or from other methods.
#' @export

round.stratEst.strategy <- function( x , digits = 0 ){
  for( i in 1:ncol(x)){
    x[,i] <- round(x[,i],digits )
  }
  return(x)
}
