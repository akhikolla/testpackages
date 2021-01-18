#' @title Compute cdf of the piecewise exponential distribution
#' @description \code{ppexp} is used to compute the cumulative distribution
#'   function of the piecewise exponential distribution. The function is ported
#'   from R to C++ via code adapted from the \code{msm} package.
#' @param q scalar. The time point at which the cdf is to be computed.
#' @param x vector or matrix. The hazard rate(s).
#' @param cuts vector. Times at which the rate changes, i.e., the interval cut
#'   points. The first element of cuts should be 0, and cuts should be in
#'   increasing order. Must be the same length as \code{x} (vector) or have the
#'   same number of columns as \code{x} (matrix)
#'
#' @details The code underlying \code{ppexp} is written in C++ and adapted from
#'   R code used in the \code{msm} package.
#'
#' @return The cumulative distribution function computed at time \code{q},
#'   hazard(s) \code{x}, and cut points \code{cuts}.
#'
#' @examples
#' # Single vector of hazard rates. Returns a single cdf value.
#' q    <- 12
#' x    <- c(0.25,0.3,0.35,0.4)
#' cuts <- c(0,6,12,18)
#' pp   <- ppexp(q,x,cuts)
#'
#' # Matrix of multiple vectors of hazard rates. Returns 10 cdf values.
#' x  <- matrix(rgamma(4*10, 0.1, 0.1), nrow=10)
#' pp <- ppexp(q,x,cuts)
#' @import methods
#' @useDynLib bayesDP
#' @export
ppexp <- function(q, x, cuts){
  if(!is.matrix(x)){
    ppout <- ppexpV(q, x, cuts)
  } else if(is.matrix(x)){
    ppout <- ppexpM(q, x, cuts)
  } else{
    stop("Error: input x is in the wrong format.")
  }
  return(ppout)
}

.onUnload <- function (libpath) {
  library.dynam.unload("bayesDP", libpath)
}
