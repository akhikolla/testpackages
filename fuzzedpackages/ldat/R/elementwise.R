
#' Apply a function to each element of an lvec
#'
#' @param x an object of type \code{\link{lvec}}.
#' @param fun the function to apply to the \code{lvec}. This function receives
#'   chunks of the lvec (which are regular R-vector) and should return a (R) 
#'   vector of the same length as its input. 
#' @param ... passed on to \code{fun}.
#'
#' @return
#' Returns a \code{link{lvec}} of the same length as the input. The type is 
#' determined by the output of \code{fun}.
#'
#' @examples
#' # Calculate square root of lvec
#' x <- as_lvec(1:10)
#' y <- elementwise(x, sqrt)
#' # of course, this is already implemented
#' sqrt(x)
#'
#' @import lvec 
#' @export
elementwise <- function(x, fun, ...) {
  chunks <- chunk(x)
  result <- NULL

  for (c in chunks) {
    d <- as_rvec(lget(x, range = c))
    r <- fun(d, ...)
    r <- as_lvec(r)
    if (is.null(result)) {
      result <- r
      length(result) <- length(x)
    } else {
      if (lvec_type(r) == 'character' && strlen(r) > strlen(result)) {
        warning("Changing maximum string length")
        strlen(result) <- strlen(r)
      }
      lset(result, range = c, r)
    }
  }
  result
}

