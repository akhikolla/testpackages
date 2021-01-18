#####
#' Simulate random variables from the Extended Empirical Saddlepoint density (ESS)
#' @description Simulate random variables from the Extended Empirical Saddlepoint density (ESS), using importance 
#'              sampling and then resampling according to the importance weights.
#'
#' @param n number of simulated vectors.
#' @param X an m by d matrix containing the data.
#' @param decay rate at which the ESS falls back on a normal density. Should be a positive number. See Fasiolo et al. (2016)
#'              for details.
#' @param ml n random variables are generated from a Gaussian importance density with covariance matrix 
#'             \code{ml*cov(X)}. By default the inflation factor is \code{ml=2}.
#' @param multicore  if TRUE the ESS densities corresponding the samples will be evaluated in parallel.
#' @param ncores   number of cores to be used.
#' @param cluster an object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster. 
#' @param ... additional arguments to be passed to \code{dsaddle}.
#' @return An n by d matrix containing the simulated vectors.
#' @details Notice that, while importance sampling is used, the output is a matrix of unweighted samples, obtained by resampling
#'          with probabilities proportional to the importance weights.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @references Fasiolo, M., Wood, S. N., Hartig, F. and Bravington, M. V. (2016). 
#'             An Extended Empirical Saddlepoint Approximation for Intractable Likelihoods. ArXiv http://arxiv.org/abs/1601.01849.
#' @examples
#' # Simulate bivariate data, where each marginal distribution is Exp(2)
#' X <- matrix(rexp(2 * 1e3), 1e3, 2)
#' 
#' # Simulate bivariate data from a saddlepoint fitted to X
#' Z <- rsaddle(1000, X, decay = 0.5)
#' 
#' # Look at first marginal distribution
#' hist( Z[ , 1] )
#' @export
#'


rsaddle <- function(n, X, decay, ml = 2,
                    multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1,  ...)
{
  prop <- rmvn(n, colMeans(X), ml*cov(X))
  
  w <- dsaddle(prop, X = X, decay = decay, multicore = multicore, ncores = ncores, cluster = cluster, ...)$llk / 
       dmvn(prop, colMeans(X), ml*cov(X))
    
  out <- prop[ sample(1:n, n, replace = TRUE, prob = w), ]
  
  return( out )
  
}