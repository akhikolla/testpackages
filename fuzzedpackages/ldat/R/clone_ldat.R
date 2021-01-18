
#' Clone an ldat object
#'
#' @param x \code{\link{ldat}} object to be cloned
#' @param ... ignored. 
#'
#' @details
#' Clones each of the vectors in the \code{\link{ldat}} object. 
#'
#' @export
clone.ldat <- function(x, ...) {
  res <- lapply(x, clone, ...) 
  attributes(res) <- attributes(x)
  res
}

