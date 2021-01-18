#' Density, cumulative density, quantiles and random number generation for the
#' extended generalized Pareto distribution 3
#' 
#' Density, cumulative density, quantiles and random number generation for the
#' EGP3 distribution of Papastathopoulos and Tawn
#' 
#' @param x,q,p Value, quantile or probability respectively.
#' @param n Number of random numbers to simulate.
#' @param kappa The power parameter (Papastathopoulos and Tawn call it the
#' shape parameter and call what we call the shape parameter the tail index.)
#' @param sigma Scale parameter.
#' @param xi Shape parameter.
#' @param u Threshold
#' @param log.d,log.p Whether or not to work on the log scale.
#' @param lower.tail Whether to return the lower tail.
#' @author Harry Southworth
#' @references I. Papastathopoulos and J. A. Tawn, Extended generalized Pareto
#' modles for tail estimation, Journal of Statistical Planning and Inference,
#' 143, 131 -- 143, 2013
#' @keywords models
#' @examples
#' 
#'   x <- regp3(1000, kappa=2, sigma=1, xi=.5)
#'   hist(x)
#'   x <- regp3(1000, kappa=2, sigma=exp(rnorm(1000, 1, .25)), xi=rnorm(1000, .5, .2))
#'   hist(x)
#'   plot(pegp3(x, kappa=2, sigma=1, xi=.5))
#' 
#' @export degp3
degp3 <- function(x, kappa=1, sigma, xi, u=0, log.d=FALSE){
    res <- log(kappa) + dgpd(x, sigma, xi, u=u, log.d=TRUE) +
        (kappa - 1) * pgpd(x, sigma, xi, u=u, lower.tail=TRUE, log.p=TRUE)

    if (log.d) res
    else exp(res)
}
