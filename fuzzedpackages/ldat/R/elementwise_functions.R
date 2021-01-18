
#' Simple elementwise functions
#'
#' These are implementations for \code{\link{lvec}} object for their regular
#' R counterparts.
#' 
#' @param x an object f type \code{\link{lvec}}
#'
#' @return 
#' Returns an \code{\link{lvec}} of the same length as the input.
#' 
#' @export
is.na.lvec <- function(x) {
  elementwise(x, is.na)
}

