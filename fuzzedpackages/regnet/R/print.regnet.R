#' print a cv.glmnet object
#'
#' Print a summary of a cv.glmnet object
#'
#' @param x cv.glmnet object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{cv.regnet}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{cv.regnet}}
#' @export
print.cv.regnet=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(x$call), "\n")
  cat("\nLambda:\n")
  print(x$lambda)
  cat("\nCV error:\n")
  print(x$mcvm)
  cat("\nPenalty:\n")
  print(x$penalty)

  # print(cbind(Df=x$df,"%Dev"=signif(x$dev.ratio,digits),Lambda=signif(x$lambda,digits)))
}
