#' Probability density function of asymmetric Laplace distributions
#'
#' \code{dald} calculates probability densities of asymmetric Laplace distributions.
#'
#' @param x Random variable.
#' @param mu Position parameter.
#' @param p Quantile.
#' @param sigma Scale parameter.
#'
#' @return probability density of \code{x}.
#' @export
#'
dald <- function(x,mu,p,sigma){
  out <- ifelse(x<=mu,p*(1-p)/sigma*exp(-(x-mu)/sigma*(p-1)),p*(1-p)/sigma*exp(-(x-mu)/sigma*p) )
  return(out)
}

#' Cumulative density function of asymmetric Laplace distributions
#'
#' \code{pald} calculates cumulative densities of asymmetric Laplace distributions.
#'
#' @param x Random variable.
#' @param mu Position parameter.
#' @param p Quantile.
#' @param sigma Scale parameter.
#'
#' @return cumulative probability density of \code{x}.
#' @export
#'
pald <- function(x,mu,p,sigma){
  out <- ifelse(x<=mu,p*exp((x-mu)*(1-p)/sigma), 1-(1-p)*exp(-(x-mu)*(p)/sigma))
  return(out)
}

#' Inverse function
#'
#' \code{inverse} generates inverse function of any given function.
#'
#' @param f pald function
#' @param mu Position parameter.
#' @param p Quantile.
#' @param sigma Scale parameter.
#' @param lower Lower bound.
#' @param upper Upper bound.
#'
#'
#' @return inversed pald
#'
#' @importFrom stats runif uniroot
#'
inverse <- function (f, mu, p, sigma, lower = -10000, upper = 10000) {
  function (y,mu,p,sigma) uniroot((function (x) f(x, mu,p,sigma) - y), lower = lower, upper = upper)[1]
}

#' Quantile function of asymmetric Laplace distributions
#'
#' \code{qald} calculates quantiles values of asymmetric Laplace distributions.
#'
#' @param y quantile value.
#' @param mu Position parameter.
#' @param p Quantile.
#' @param sigma Scale parameter.
#'
#' @return quantile value.
#' @export
#'
qald <- inverse(pald)

#' Random number generator of asymmetric Laplace distributions
#'
#' \code{rald} generates random numbers from asymmetric Laplace distributions.
#'
#' @param n Number of random numbers to be generated.
#' @param mu Position parameter.
#' @param p Quantile.
#' @param sigma Scale parameter.
#'
#' @return random numbers.
#' @export
#'
rald <- function(n,mu,p,sigma){
  tmp <- runif(n)
  out <-  NA
  for (i in 1:n){
    out[i] <- qald(tmp[i],mu,p,sigma)
  }
  return(unlist(out))
}
