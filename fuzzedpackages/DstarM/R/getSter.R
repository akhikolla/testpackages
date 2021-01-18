#' Estimate variance of nondecision density
#'
#' @param res An object of class \code{D*M}.
#'
#' @details The object \code{res} can either be output from \code{estDstarM} or output from \code{estND}.
#' if the former is supplied, \code{getSter} attempts to calculate the variance of the
#' nondecision distribution by subtracting the variance of the model distribution from the
#' variance of the data distribution. If the latter is supplied, the variance is calculated by
#' integrating the nondecision distribution.
#'
#'

#' @export
# works for nondecision densities, check if works for estDstarM output.
getSter <- function(res) {
  
  if (is.DstarM.fitND(res)) {
    return(apply(res$r.hat, 2, nth.cmomentS, x = res$tt, nth = 2))
  } else if (is.DstarM.fitD(res)) {
    return(res$var.dat - res$var.m)
  } else {
    stop("res should be output from either estDstarM or estND.")
  }
}
