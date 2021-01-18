
#' Read elements from an ldat object
#'
#' @param x the \code{\link{ldat}} to read from
#' @param ... passed on to \code{\link{lget.lvec}}.
#'
#' @details
#' Indexing using \code{index} should follow the same rules as indexing a regular
#' data.frame using a logical or numeric index. The range given by \code{range} 
#' includes both end elements. So, a range of \code{c(1,3)} selects the first 
#' three elements. 
#'
#' @return
#' Returns an \code{\link{ldat}} with the selected elements. In order to convert
#' the selection to an R-vector \code{\link{as_rvec}} can e used.
#'
#' @export
lget.ldat <- function(x, ...) {
  res <- lapply(x, lget, ...)
  structure(res, class = 'ldat')
}

