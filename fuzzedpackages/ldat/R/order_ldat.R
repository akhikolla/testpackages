
#' Order an ldat
#'
#' @param x \code{\link{ldat}} to sort
#' @param ... unused.
#'
#' @return
#' Returns the order of \code{x}. Unlike the default \code{\link{order}} 
#' function in R, the sort used is not stable (e.g. in case there are multiple
#' records with the same value in \code{x}, there relative order after sorting
#' is not defined). 
#'
#' @examples
#' x <- as_ldat(iris)
#' o <- order(x[c("Sepal.Width", "Sepal.Length")])
#'
#' @useDynLib ldat
#' @export
order.ldat <- function(x, ...) {
  if (!is_ldat(x))
    stop("x should be of type ldat")
  #o <- .Call("order_ldat", x, PACKAGE = "ldat")
  o <- order_ldat_cpp(x)
  structure(o, class = "lvec")
}


