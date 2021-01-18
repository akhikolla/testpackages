#' print a springer result
#'
#' Print a springer result
#'
#' @param x springer result
#' @param digits significant digits in printout.
#' @param ... other print arguments
#' @seealso \code{\link{springer}}
#' @export
print.springer=function(x, digits = max(3, getOption("digits") - 3),...){
  cat("\nCoefficients:\n")
  print(x$coefficient, digits)
}

