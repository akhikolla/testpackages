#' Title
#'
#' @param x
#' @param y
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
KdCov = function(x,y, sigma){
  n <- nrow(as.matrix(x))
    if (length(y)!=n)
        stop( "x and y must be the same size")
    x=sqrt(2-2*exp(-as.matrix(dist(as.matrix(x),diag=TRUE, upper=TRUE))/sigma))
    y <- as.matrix(dist(as.matrix(y),diag=TRUE, upper=TRUE))
    y[which(y!=0)]=1

    return(dcov(x,y))
}
