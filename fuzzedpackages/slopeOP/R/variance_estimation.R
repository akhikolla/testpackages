
######################################################
#############       sdHallDiff       #################
######################################################

#' sdHallDiff
#' @description Estimation of the standard deviation using the HallDiff estimator
#' @param data vector of data to segment: a univariate time series
#' @return an estimation of the sd
#' @examples
#' myData <- slopeData(index = c(1,100,200,300), states = c(0,5,3,6), noise = 1)
#' sdHallDiff(data = myData)

sdHallDiff <- function(data)
{
  wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
  corrector <- wei[4]^2 + (wei[3]-wei[4])^2 + (wei[2]-wei[3])^2 + (wei[1]-wei[2])^2 + wei[1]^2

  z <- diff(data) #diff data
  n <- length(z)
  mat <- wei %*% t(z)
  mat[2, -n] <- mat[2, -1]
  mat[3, -c(n-1, n)] <- mat[3, -c(1, 2)]
  mat[4, -c(n-2, n-1, n)] <- mat[4, -c(1, 2, 3)]
  sd <- sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2)/((n-3)*corrector))

  return(sd)
}
