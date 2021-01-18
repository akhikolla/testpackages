#' Density, cumulative density, quantiles and random number generation for the
#' generalized extreme value distribution
#' 
#' Density, cumulative density, quantiles and random number generation for the
#' generalized extreme value distribution
#' 
#' Random number generation is done as a transformation of the Gumbel
#' distribution; Gumbel random variates are generated as the negative logarithm
#' of standard exponentials.
#' 
#' @param x,q,p Value, quantile or probability respectively.
#' @param n Number of random numbers to simulate.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#' @param xi Shape parameter.
#' @param log.d,log.p Whether or not to work on the log scale.
#' @param lower.tail Whether to return the lower tail.
#' @author Harry Southworth
#' @keywords models
#' @examples
#' 
#'   x <- rgev(1000, mu=0, sigma=1, xi=.5)
#'   hist(x)
#'   x <- rgev(1000, mu=0, sigma=exp(rnorm(1000, 1, .25)), xi=rnorm(1000, .5, .2))
#'   hist(x)
#'   plot(pgev(x, mu=0, sigma=1, xi=.5))
#' 
#' @export dgev
dgev <- function(x, mu, sigma, xi, log.d=FALSE){
  ## shift and scale
  x <- (x - mu) / sigma

  xix <- .specfun.safe.product(xi, x)
  logrel <- .log1prel(xix) * x

  log.density <- -log(sigma) - log1p(xix) - logrel - exp(-logrel)
  ## make exp(Inf) > Inf
  log.density[logrel==(-Inf)] <- -Inf

  if (!log.d) {
    exp(log.density)
  } else {
    log.density
  }
}
