
#' Sort an ldat
#'
#' @param x \code{\link{ldat}} to sort
#' @param decreasing unused (a value unequal to \code{FALSE} will generate an error).
#' @param ... unused.
#'
#' @return
#' Sorts \code{x} and returns a sorted copy of \code{x}. 
#'
#' @examples
#' x <- as_ldat(iris)
#' sort(x)
#'
#' @export
sort.ldat <- function(x, decreasing = FALSE, ...) {
  if (decreasing != FALSE) stop("decreasing is not supported yet.")
  if (!is_ldat(x)) stop("x should be of type ldat.")
  o <- order(x)
  x <- x[o, ]
  x
}

