#' Random generator from the beta distribution
#'
#' The function randomly generates values from the beta distribution with a mean-precision parameterization.
#' @param n the number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu the mean parameter of the beta distribution. It must lie in (0, 1).
#' @param phi the precision parameter of the Beta distribution. It must be a positive real value.
#'
#' @return A vector of length \code{n}.
#'
#' @examples rBeta_mu(100, mu = 0.5, phi = 30)
#'
#' @references{
#' Ferrari, S.L.P., and Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. Journal of Applied Statistics, \bold{31}(7), 799--815. doi:10.1080/0266476042000214501
#' }
#'
#' @import stats
#'
#' @export
#'

rBeta_mu <- function(n, mu, phi){
  if (length(n)>1) n <- length(n)
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  alpha1 <- mu*phi
  alpha2 <- (1-mu)*phi
  n <- floor(n)
  return(rbeta(n, shape1=alpha1, shape2=alpha2))
}

#' Random generator from the flexible beta distribution
#'
#' The function randomly generates values from the flexible beta distribution.
#' @param n the number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu the mean parameter of the flexible beta distribution. It must lie in (0, 1).
#' @param phi the precision parameter of the Flexible Beta distribution. It must be a positive real value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param w the normalized distance among clusters. It must lie in (0, 1).
#'
#' @return A vector of length  \code{n}.
#'
#' @examples rFB(100,0.5,30,0.3,0.6)
#'
#' @references {
#' Migliorati, S., Di Brisco, A. M., Ongaro, A. (2018). A New Regression Model for Bounded Responses. Bayesian Analysis, \bold{13}(3), 845--872. doi:10.1214/17-BA1079
#' }
#'
#' @import stats
#'
#' @export

rFB <- function(n, mu, phi, p, w){
  if (length(n)>1) n <- length(n)
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(w < 0 | w > 1)) stop("Parameter w has to be between 0 and 1")
  n <- floor(n)
  wtilde <- w*min(mu/p, (1-mu)/(1-p))
  lambda1 <- mu + (1-p)*wtilde
  lambda2 <- mu-p*wtilde
  v <- rbinom(n,1,prob=p)
  x <- vector(mode="numeric", length = n)
  x[v==1] <- rBeta_mu(length(which(v==1)),lambda1,phi)
  x[v==0] <- rBeta_mu(length(which(v==0)),lambda2,phi)
  return(x)
}



#' Random generation from the variance-inflated beta distribution
#'
#' The function randomly generates values from the variance-inflated beta distribution.
#' @param n the number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mu the mean parameter of the Variance-Inflated distribution. It must lie in (0, 1).
#' @param phi the precision parameter of the Variance-Inflated distribution. It must be a positive real value.
#' @param p the mixing weight. It must lie in (0, 1).
#' @param k the extent of the variance inflation. It must lie in (0, 1).
#'
#' @return A vector of length  \code{n}.
#'
#' @examples rVIB(100,0.5,30,0.3,0.6)
#'
#' @references{
#'  Di Brisco, A. M., Migliorati, S., Ongaro, A. (2020) Robustness against outliers: A new variance inflated regression model for proportions. Statistical Modelling, \bold{20}(3), 274--309. doi:10.1177/1471082X18821213
#' }
#'
#' @import stats
#'
#' @export

rVIB <- function(n, mu, phi, p, k){
  if (length(n)>1) n <- length(n)
  if (any(mu < 0 | mu > 1)) stop("Parameter mu has to be between 0 and 1")
  if (any(phi < 0)) stop("Parameter phi has to be greater than 0")
  if (any(p < 0 | p > 1)) stop("Parameter p has to be between 0 and 1")
  if (any(k < 0 | k > 1)) stop("Parameter k has to be between 0 and 1")
  n <- floor(n)
  v <- rbinom(n,1,prob=p)
  x <- vector(mode="numeric", length = n)
  x[v==1] <- rBeta_mu(length(which(v==1)),mu,phi*k)
  x[v==0] <- rBeta_mu(length(which(v==0)),mu,phi)

  return(x)
}

