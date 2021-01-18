#' Title
#'
#' @param x
#' @param y
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
dCor <- function(x, y, alpha) {
  if (missing(alpha))
    alpha <-1
  n <- nrow(as.matrix(x))
  if (length(y)!=n)
    stop( "x and y must be the same size")
  x <- as.matrix(dist(as.matrix(x),diag=TRUE, upper=TRUE))
  y <- as.matrix(dist(as.matrix(y),diag=TRUE, upper=TRUE))
  y[which(y!=0)]=1
  return(dcor(x,y,alpha))
}
