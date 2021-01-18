#####
#' Evaluate the density of a multivariate Gaussian fit
#' @description Given a sample X, it gives a pointwise evaluation of the multivariate normal (MVN) density fit at position y.
#'
#' @param y points at which the MVN is evaluated. It can be either a d-dimensional vector or an n by d matrix, 
#'        each row indicating a different position.
#' @param X an n by d matrix containing the data.
#' @param log if TRUE the log-density is returned.
#' @param verbose currently not used.
#' @param alpha tuning parameter of \code{robCov}, see \code{?robCov} for details.
#' @param beta tuning parameter of \code{robCov}, see \code{?robCov} for details.
#' @return A vector where the i-th entry is the density corresponding to the i-th row of y.
#' @details The covariance matrix is estimated robustly, using the \code{robCov} function.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon N. Wood.
#' @examples 
#' library(esaddle)
#' X <- matrix(rnorm(2 * 1e3), 1e3, 2) # Sample used to fit a multivariate Gaussian
#' demvn(rnorm(2), X, log = TRUE)      # Evaluate the fitted log-density at a random location
#' @export demvn
#'
demvn <- function(y, X, log = FALSE, verbose = TRUE, alpha=2, beta=1.25)
{
  
  if( !is.matrix(y) ) y <- matrix(y, 1, length(y))
  
  tmp <- robCov( t(X), alpha = alpha, beta = beta )
  
  # If there are some statistics with zero variace we remove them
  if( length(tmp$lowVar) ) stop("The columns of X indexed ", tmp$lowVar, " have zero variance.") #y <- y[-tmp$lowVar]
  
  llk <- apply(y, 
               1, 
               function(input) .demvn(y = input, L = tmp, log = log) )
  
  return(llk)
}


.demvn <- function(y, L, log)
{
  rss <- sum( (L$E%*%as.vector(y-L$mY))^2 )
  
  llk <- -rss/2 - L$half.ldet.V - log(2 * pi) * length(y) / 2
  
  if( !log ) llk <- exp( llk )
  
  return(llk)
}