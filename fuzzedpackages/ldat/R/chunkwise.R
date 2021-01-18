
#' Process an lvec in chunks
#'
#' @param x the \code{\link{lvec}}.
#' @param init initialisation function. This function should accept an 
#'   \code{\link{lvec}} as its first argument and return an initial value
#'   for the state.
#' @param update update function. Called for each chunk of data. Receives 
#'   the current value of the state as its first argument and the next chunk
#'   of data as its second argument. Should return an updated state. This 
#'   function can be called multiple times.
#' @param final finaliser function. Is called after processing the complete
#'   lvec. Receives the final state as its first argument. Should return the 
#'   end result.
#' @param ... optional arguments passed on to the supplied functions.
#'
#' @details
#' For examples of its use see \code{\link{mean.lvec}} and 
#' \code{\link{sum.lvec}}.
#'
#' @import lvec
#' @export
chunkwise <- function(x, init, update, final, ...) {
  chunks <- chunk(x)
  state <- init(x, ...)
  for (c in chunks) {
    d <- as_rvec(lget(x, range = c))
    state <- update(state, d, ...)
  }
  final(state, ...)
}

