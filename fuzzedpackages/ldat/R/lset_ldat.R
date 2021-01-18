
#' Set values is an ldat object
#'
#' @param x an object of type \code{\link{ldat}}
#' @param index a numeric of logical vector with indices at which 
#'   the values should be set. 
#' @param values a vector with new values.
#' @param range a numeric vector of length 2 specifying a range of elements 
#'   to select. Specify either \code{index} or \code{range}. 
#' @param ... ignored.
#'
#' @details
#' When values is a vector the values are assigned to each column in
#' \code{x}. Otherwise, vectors is assumed to be a list or data.frame
#' of the same length as \code{x}. Each element of \code{values} is 
#' assigned to the corresponding element of \code{x}.
#'
#' @export
lset.ldat <- function(x, index = NULL, values, range = NULL, ...) {
  if (!is.null(ncol(values))) {
    # data.frame like object
    stopifnot(ncol(values) > 0)
    j <- 1
    for (i in seq_len(ncol(x))) {
      x[[i]] <- lset(x[[i]], values[[j]], index = index, range = range)
      j <- (j %% ncol(values)) + 1
    }
  } else {
    # assume vector: assign to each column
    lapply(x, lset, index = index, range = range, values = values)
  }
  x
}

