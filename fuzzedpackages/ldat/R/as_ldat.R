
## TODO: tests

#' Convert r-objects to \code{\link{ldat}}'s
#' 
#' @param x object to convert
#' @param ... further arguments passed to or from other methods
#'
#' @return
#' Returns a \code{ldat} with columns corresponding to the columns in \code{x}.
#' When \code{x} is not a \code{data.frame} it is first converted to a 
#' \code{data.frame} using a call to \code{\link{as.data.frame}}.
#'
#' @examples
#' a <- data.frame(a = 1:10, b = rnorm(10))
#' b <- as_ldat(a)
#'
#' @export
as_ldat <- function(x, ...) {
  UseMethod("as_ldat")
}

#' @rdname as_ldat
#' @export
as_ldat.data.frame <- function(x, ...) {
  do.call(ldat, x)
}

#' @rdname as_ldat
#' @export
as_ldat.default <- function(x, ...) {
  x <- as.data.frame(x)
  as_ldat(x, ...)
}

#' @rdname as_ldat
#' @export
as_ldat.ldat <- function(x, ...) {
  x
}
