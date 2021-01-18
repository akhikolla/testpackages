#' Density function
#'
#' @param rt vector of reaction times
#' @param tt grid to evaluate the density on
#'
#' @return a vector of \code{length(tt)}
#'
#' @details Can be passed to the argument \code{densityMethod} of \code{\link{estDstarM}}. This function is a minimal
#' example to use as custom smoothing function.
#'
#' @examples
#' x <- rgamma(1e5, 1, 1)
#' tt <- seq(0, 5, .01)
#' d <- Density(x, tt)
#' hist(x, freq = FALSE)
#' lines(tt, DstarM:::Density(x, tt))

#' @export
Density <- function(rt, tt) {
  return(stats::density(x = rt, from = min(tt), to = max(tt), n = length(tt))$y)
}
