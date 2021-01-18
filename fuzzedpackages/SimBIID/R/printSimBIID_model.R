#' @title Prints \code{SimBIID_model} objects
#'
#' @description Print method for \code{SimBIID_model} objects.
#'
#' @param x    A \code{SimBIID_model} object.
#' @param ...           Not used here.
#'
#' @return Prints parsed \code{Rcpp} code to the screen.
#' 
#' @export

print.SimBIID_model <- function(x, ...) {
    writeLines(x$code)
}