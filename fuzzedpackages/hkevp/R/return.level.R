#' @title The associated return level
#' @author Quentin Sebille
#'
#'
#' @description 
#' Computation of the associated return level with given period and GEV parameters.
#' 
#' 
#' 
#' @param period
#' An integer indicating the wished return period T.
#' 
#' @param loc
#' A numerical value or vector for the GEV location parameter. Must be of length one or same length as \code{scale} and/or \code{shape}.
#' 
#' @param scale
#' A numerical value or vector for the GEV scale parameter. Must be of length one or same length as \code{loc} and/or \code{shape}.
#' 
#' @param shape
#' A numerical value or vector for the GEV shape parameter. Must be of length one or same length as \code{loc} and/or \code{scale}.
#'
#'
#'
#' @details
#' The \eqn{T}-year return level is a common value of risk in Extreme Value Theory. It represents the value that is expected to be exceeded once over \eqn{T} years by the annual maxima. Given the parameters \eqn{\mu}, \eqn{\sigma} and \eqn{\xi} of the GEV distribution associated to the yearly maxima, we can compute the associated \eqn{T}-return level \eqn{y_T} by:
#' \deqn{y_T := \mu + \frac{\sigma}{\xi} \left[ \log\left(\frac{T}{T-1}\right)^{-\xi} -1 \right] ~.}
#'
#'
#'
#' @return 
#' A numerical value or a numerical vector, depending on the input arguments \code{loc}, \code{scale}, \code{shape}
#'
#'
#' @examples
#' return.level(period = 100, loc = 1, scale = 1, shape = 1)
#' return.level(period = 200, loc = 1:10, scale = 1, shape = 0)
#' 
#' 
return.level <- function(period, loc, scale, shape) {
  ## Errors
  if ((period != round(period)) || (period <= 1)) stop("'period' must be an integer greater than 1!")
  LENGTHS <- c(length(loc), length(scale), length(shape))
  if (!all(LENGTHS==1 | LENGTHS==max(LENGTHS))) stop("GEV parameters must have same length (or length one)!")
  
  ## Initialisation of the result table
  loc <- as.vector(loc)
  scale <- as.vector(scale)
  shape <- as.vector(shape)
  gev <- cbind(loc, scale, shape)
  result <- rep(NA, nrow(gev))
  
  ## Separation of two cases : shape == 0 and shape != 0
  isshape0 <- gev[,3] == 0
  result[isshape0] <- gev[isshape0,1] - gev[isshape0,2]*log(log(period/(period - 1)))
  result[!isshape0] <- gev[!isshape0,1] + gev[!isshape0,2]*(log(period/(period - 1)) ^ (-gev[!isshape0,3]) - 1)/gev[!isshape0,3]
  
  return(result)
}
