#' Azzalini skew logistic distribution
#'
#' Generates random numbers from the skew logistic distribution
#'
#' @inheritParams sn::rsc
#' @return Vector of random numbers from Azzalini skew logistic distribution.
#' @export
rsl <- function(n=1, xi=0, omega=1, alpha=0, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(alpha))
      stop("You cannot set both 'dp' and component parameters")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
  }
  u1 <- rlogis(n)
  u2 <- rlogis(n)
  id <- (u2 > alpha*u1)
  u1[id] <- (-u1[id])
  z <- u1
  y <- xi+omega*z
  attr(y, "family") <- "SL"
  attr(y, "parameters") <- c(xi,omega,alpha)
  return(y)
}


#' Mixture of 2 normal distributions
#'
#' Generates random numbers from a mixture of 2 normal distributions
#'
#' @param n number of observations
#' @param mean1 mean of the first normal
#' @param sd1 standard deviation of the first normal
#' @param mean2 mean of the second normal
#' @param sd2 standard deviation of the second normal
#' @param p probability of the first normal
#' @return Vector of random numbers from the specified mixture of normals.
#' @export
rmixnorm <- function(n, mean1 = 0, sd1 = 1, mean2 = 0, sd2 = 1, p = 0.5)
{
  take_first <- sample(0:1, n, replace = TRUE, prob = c(1-p, p))
  take_first*rnorm(n, mean1, sd1) + (!take_first)*rnorm(n, mean2, sd2)
}
