
#' Select a range of records from an object
#' 
#' @param x the object to select items from
#' @param range a numeric vector with two elements specifying the range to 
#'   select. 
#' @param begin the first element to select.
#' @param end the last element to select.
#' @param as_r convert the result to an R-object. 
#' @param ... ignored; used to pass additional arguments to other methods. 
#' 
#' @examples
#' x <- as_lvec(1:20)
#' # Select elements 5:7
#' slice_range(x, range = c(5, 7))
#' slice_range(x, begin = 5, end = 7)
#' slice_range(x, range = c(5, 10), end = 7)
#' # also works for R-vectors
#' slice_range(1:20, range = c(5,7))
#' # convert lvec to rvec
#' slice_range(x, range = c(5,7), as_r = TRUE)
#' 
#' @rdname slice_range
#' @export
slice_range <- function(x, range, begin = range[1], end = range[2], ...) {
  UseMethod("slice_range")
}


#' @rdname slice_range
#' @export
slice_range.lvec <- function(x, range, begin = range[1], end = range[2], 
    as_r = FALSE, ...) {
  res <- lget(x, range = c(begin, end))
  if (as_r) as_rvec(res) else res
}

#' @rdname slice_range
#' @export
slice_range.ldat <- function(x, range, begin = range[1], end = range[2],
    as_r = FALSE, ...) {
  res <- lapply(x, slice_range, begin = begin, end = end, as_r = as_r, ...)
  if (as_r) data.frame(res, stringsAsFactors = FALSE) else structure(res, class = 'ldat')
}

#' @rdname slice_range
#' @export
slice_range.default <- function(x, range, begin = range[1], end = range[2], 
    ...) {
  if (begin < 1 || end > length(x)) stop('Index out of range')
  x[seq.int(begin, end, by = 1)]
}


#' @rdname slice_range
#' @export
slice_range.data.frame <- function(x, range, begin = range[1], end = range[2], 
    ...) {
  if (begin < 1 || end > nrow(x)) stop('Index out of range')
  x[seq.int(begin, end, by = 1), , drop = FALSE]
}

