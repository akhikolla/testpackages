#' Sample from seemingly unrelated regression
#'
#' This function is a wrapper function that performs either (1) Direct Monte Carlo or (2) Gibbs sampling
#' of the SUR model depending on whether 1 or 2 data sets are specified.
#' 
#' @include sur_sample_dmc.R
#' @include sur_sample_powerprior.R
#' @importFrom stats as.formula lm model.matrix
#' 
#' @param formula.list A list of formulas, each element giving the formula for the corresponding endpoint.
#' @param data A \code{data.frame} containing all the variables contained in \code{formula.list]}
#' @param M Number of samples to be drawn
#' @param histdata A \code{data.frame} of historical data to apply power prior on
#' @param Sigma0 \emph{optional} a \eqn{J \times J} \code{matrix} giving the initial covariance matrix. Default is the MLE. Ignored if \code{histdata} is null
#' @param a0 A scalar between 0 and 1 giving the power prior parameter. Ignored if \code{histdata} is null
#' @param burnin A non-negative integer giving the burn-in parameter. Ignored if \code{histdata} is null as direct Monte Carlo is performed
#' @param thin A positive integer giving the thin parameter. Ignored if \code{histdata} is null as direct Monte Carlo is performed
#' @return A list. First element is posterior draws. Second element is list of JxJ covariance matrices.
#' 
#' @examples
#' ## Taken from bayesm package
#' if(nchar(Sys.getenv("LONG_TEST")) != 0) {M=1000} else {M=10}
#' set.seed(66)
#' ## simulate data from SUR
#' beta1 = c(1,2)
#' beta2 = c(1,-1,-2)
#' nobs = 100
#' nreg = 2
#' iota = c(rep(1, nobs))
#' X1 = cbind(iota, runif(nobs))
#' X2 = cbind(iota, runif(nobs), runif(nobs))
#' Sigma = matrix(c(0.5, 0.2, 0.2, 0.5), ncol = 2)
#' U = chol(Sigma)
#' E = matrix( rnorm( 2 * nobs ), ncol = 2) %*% U
#' y1 = X1 %*% beta1 + E[,1]
#' y2 = X2 %*% beta2 + E[,2]
#' X1 = X1[, -1]
#' X2 = X2[, -1]
#' data = data.frame(y1, y2, X1, X2)
#' names(data) = c( paste0( 'y', 1:2 ), paste0('x', 1:(ncol(data) - 2) ))
#' ## run DMC sampler
#' formula.list = list(y1 ~ x1, y2 ~ x2 + x3)
#' 
#' ## Fit models
#' out_dmc = sur_sample( formula.list, data, M = M )            ## DMC used
#' out_powerprior = sur_sample( formula.list, data, M, data )   ## Gibbs used
#' @export
sur_sample <- function(
  formula.list, 
  data, 
  M, 
  histdata = NULL,
  Sigma0 = NULL, 
  a0 = 1, 
  burnin = 0, 
  thin = 1
) {
  
  if ( is.null( histdata ) ) {
    message('Direct Monte Carlo sampling used')
    return ( sur_sample_dmc( formula.list, data, M ) )
  } else {
    message('Gibbs sampling used for power prior model')
    return ( sur_sample_powerprior( formula.list, data, histdata, M, Sigma0, a0, burnin, thin ) )
  }
}

