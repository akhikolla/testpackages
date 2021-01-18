#' Plot outliers detected by robets model
#' 
#' @param object An object of class \code{robets}.
#' @param xlab Label of the x-axis.
#' @param ylab Label of the y-asix.
#' @param type Character indicating the type of plot, just as in \code{plot}.
#' @param ... Other plotting parameters.
#'
#' @examples
#' model <- robets(nottem)
#' plotOutliers(model)
#' @seealso \code{\link{plot.robets}}
#' @export
plotOutliers <- function(object, xlab="", ylab="", type="l", ...) {
  y <- object$x
  ny <- length(y)
  y2 <- na.contiguous(y)
  if(ny != length(y2))
    warning("Missing values encountered. Using longest contiguous portion of time series, just as in robets.")
  y <- as.ts(y2)
  n <- length(y)
  plot(y, xlab = xlab, ylab = ylab, type = type, main=paste("Outlier detection with", object$method,"method"), ...)
  
  m <- frequency(y)
  srt <- start(y)
  
  for( t in 1:n){
    if(object$outlier[t]){
      points(srt[1]+(t-1+srt[2]-1)/m, y[t],col="red",pch = 19)
    }
  }
  
  if(sum(object$outlier) == 0) message("No outliers detected.")
}