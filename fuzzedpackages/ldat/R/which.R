
#' @rdname which
#' @export
which <- function(x, ...) {
  UseMethod("which")
}

#' @rdname which
#' @export
which.default <- function(x, ...) {
  base::which(x, ...)
}

#' Give the `TRUE` indices of an lvec
#'
#' @param x logical \code{\link{lvec}} to get the indices from
#' @param ... not used
#' 
#' @return
#' Returns a numeric lvec with the indices of the elements of \code{x} that are 
#' TRUE. 
#'
#' @examples
#' x <- as_lvec(runif(1E6) > 0.1)
#' which(x)
#'
#' @rdname which
#' @export
which.lvec <- function(x, ...) {
  chunks <- chunk(x)
  result <- NULL
  for (c in chunks) {
    d <- which(as_rvec(x[range = c]))
    d <- d + (c[1] - 1)
    result <- append.lvec(result, d, clone = FALSE)
  }
  result
}
