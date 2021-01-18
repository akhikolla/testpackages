#' print a roben object
#'
#' Print a summary of a roben object
#'
#' @param x roben object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{roben}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{roben}}
#' @export
print.roben=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(x$call), "\n")
  # cat("\nLambda:\n")
  # print(x$lambda)
  cat("\nCoefficients:\n")
  print(x$coefficient, digits)
  cat("Class:\n")
  print(class(x))
}


#' print a GxESelection object
#'
#' Print a summary of a GxESelection object
#'
#' @param x GxESelection object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{GxESelection}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{GxESelection}}
#' @export
print.GxESelection=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nMethod:\n")
  print(x$method)
  cat("\n")
  print(x$summary)
}


#' print a roben.pred object
#'
#' Print a summary of a roben.pred object
#'
#' @param x roben.pred object.
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @usage \method{print}{roben.pred}(x, digits = max(3, getOption("digits") - 3), \dots)
#' @seealso \code{\link{predict.roben}}
#' @export
print.roben.pred=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nPMSE:\n")
  print(x$error, digits)
  cat("\npredicted ", length(x$y.pred), " y (list component y.pred)", sep = "")
}
