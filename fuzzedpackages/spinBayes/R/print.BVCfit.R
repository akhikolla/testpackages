#' print a BVCfit object
#'
#' Print a summary of a BVCfit object
#'
#' @param x BVCfit object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{BVCfit}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{BVCfit}}
#' @export
print.BVCfit=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(x$call), "\n")
  # cat("\nLambda:\n")
  # print(x$lambda)
  cat("\nCoefficients:\n")
  print(x$coefficient)
  cat("Class:\n")
  print(class(x))

  # print(cbind(Df=x$df,"%Dev"=signif(x$dev.ratio,digits),Lambda=signif(x$lambda,digits)))
}


#' print a BVCfit.pred object
#'
#' Print a summary of a BVCfit.pred object
#'
#' @param x BVCfit object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{BVCfit.pred}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{predict.BVCfit}}
#' @export
print.BVCfit.pred=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nPMSE:\n")
  print(x$pmse)
  cat("\npredicted ", length(x$y.pred), " y (list component y.pred)", sep = "")
}


#' print a BVSelection object
#'
#' Print a summary of a BVSelection object
#'
#' @param x BVSelection object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{BVSelection}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{BVSelection}}
#' @export
print.BVSelection=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nMethod:\n")
  print(x$method)
  cat("\n")
  print(x$summary)
}
