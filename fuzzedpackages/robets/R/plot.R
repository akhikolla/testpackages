#' Plot robets model
#' 
#' @param x An object of class \code{robets}.
#' @param ... Other plotting parameters.
#' @method plot robets
#'
#' @examples
#' model <- robets(nottem)
#' plot(model)
#' @seealso \code{\link{plotOutliers}, \link{plot.ets}}
#' @export
plot.robets <- function(x,...)
{
  if(!is.null(x$lambda))
    y <- BoxCox(x$x,x$lambda)
  else
    y <- x$x
  if(x$components[3]=="N" & x$components[2]=="N")
  {
    plot(cbind(observed=y, level=x$states[,2]),
         main=paste("Decomposition by",x$method,"method"),...)
  }
  else if(x$components[3]=="N")
  {
    plot(cbind(observed=y, level=x$states[,2], slope=x$states[,"b"]),
         main=paste("Decomposition by",x$method,"method"),...)
  }
  else if(x$components[2]=="N")
  {
    plot(cbind(observed=y, level=x$states[,2], season=x$states[,"s1"]),
         main=paste("Decomposition by",x$method,"method"),...)
  }
  else
  {
    plot(cbind(observed=y, level=x$states[,2], slope=x$states[,"b"],
               season=x$states[,"s1"]),
         main=paste("Decomposition by",x$method,"method"),...)
  }
}