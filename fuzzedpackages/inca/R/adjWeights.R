#' @title
#' Function for Weights Adjustments
#'
#' @description
#' This function provides a trimming procedure to force the weights to be within the provided boundaries
#' 
#' @param weights A numerical vector of weights
#' @param lower A numerical vector of lower bounds
#' @param upper A numerical vector of upper bounds
#' 
#' @details
#' The function produces trimmed weights, which will be the input for the rounding
#' technique before integer calibration. When the weights are bounded, the function rounds-up 
#' the lower bounds and rounds-down the upper. If the condition \code{upper > lower + 1},
#' an error is returned.
#' 
#' @return 
#' A vector of adjusted weights
#' 
#' @examples
#' library(inca)
#' w <- rnorm(150, 0, 2)
#' aw <- adjWeights(w, runif(150, -3, -1), runif(150, 1, 3))
#' hist(aw, main = "Adjusted weights")
#' 
#' @export

"adjWeights" <- function(weights, lower = -Inf, upper = +Inf) {
  mylower <- ceiling(lower)
  myupper <- floor(upper)
  if (any(mylower > myupper)) stop("lower bounds cannot be greater than the upper bounds.")
  res <- pmin(pmax(weights, mylower), myupper)
  return(res)
}
