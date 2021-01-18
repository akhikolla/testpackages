#' Summary robets model
#' 
#' @param object An object of class \code{robets}.
#' @param ... Other undocumented arguments.
#' 
#' @return A number of training set error measures: ME (mean error), RMSE (root mean squared error), MAE (mean absolute error), MPE (mean percentage error), MAPE (mean absolute percentage error), MedianE (median error), RTSE (root tau squared error), RTSPE (root tau squared percentage error).
#'
#' @examples
#' model <- robets(nottem)
#' summary(model)
#' @export
summary.robets <- function(object, ...)
{
  res <- object$x - object$fitted # one step ahead prediction errors, in sample
  pe <- res/object$fitted * 100
  me <- mean(res,na.rm=TRUE)
  mse <- mean(res^2,na.rm=TRUE)
  mae <- mean(abs(res),na.rm=TRUE)
  mape <- mean(abs(pe),na.rm=TRUE)
  mpe <-  mean(pe,na.rm=TRUE)
  errors <- c(me, sqrt(mse), mae, mpe, mape)
  names(errors) <- c("ME", "RMSE", "MAE", "MPE", "MAPE")
  
  print(object)
  cat("Training set: non robust error measures: \n")
  print(errors)
  
  cat("Training set: robust error measures: \n")
  me <- median(res,na.rm=TRUE)
  t2 <- tau2(res)
  t2p <- tau2(pe)
  errorsr <- c(me, t2, t2p)
  names(errorsr) <- c("MedianE", "RTSE", "RTSPE")
  print(errorsr)
  
  return(c(errors,errorsr))
}