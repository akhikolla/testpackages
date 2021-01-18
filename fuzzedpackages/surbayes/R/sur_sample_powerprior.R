#' Sample from SUR posterior via power prior
#'
#' This function uses Gibbs sampling to sample from the posterior density of a SUR model
#' using the power prior.
#' 
#' @importFrom stats as.formula lm model.matrix terms
#' 
#' @param formula.list A list of formulas, each element giving the formula for the corresponding endpoint.
#' @param data A \code{data.frame} containing all the variables contained in \code{formula.list]}
#' @param histdata A \code{data.frame} of historical data to apply power prior on
#' @param M Number of samples to be drawn
#' @param Sigma0 A \eqn{J \times J} \code{matrix} giving the initial covariance matrix. Default is the MLE.
#' @param a0 A scalar between 0 and 1 giving the power prior parameter
#' @param burnin A non-negative integer giving the burn-in parameter
#' @param thin A positive integer giving the thin parameter
#' 
#' @return A list. First element is posterior draws. Second element is list of JxJ covariance matrices.
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
#' ## fit using historical data as current data set--never done in practice
#' out = sur_sample_powerprior( formula.list, data, histdata = data, M = M )
#' @export
sur_sample_powerprior <- function(
  formula.list, data, histdata, M, Sigma0 = NULL, a0 = 1, burnin = 0, thin = 1
) {
  
  ##
  ## Perform checks
  ## 
    ## check formula.list
    if ( any( sapply( formula.list, class ) != 'formula' ) ) {
      stop('formula.list must be a list of formulas')
    }
    
    ## check data and histdata
    if ( class(data) != 'data.frame' ) { stop('data must be a data.frame') }
    if ( class(histdata) != 'data.frame' ) { stop('histdata must be a data.frame') }
  
    ## Check if M is positive integer
    if ( M <= 0 | M %% 1 != 0 ) { stop('M must be a positive integer') }
  
    ## Make sure Sigma0 is NULL or a matrix of proper dimensions
    if ( !is.null(Sigma0) ) {
      if ( class(Sigma0) != 'matrix' ) {
        stop('Sigma0 must be a matrix')
      }
      if ( nrow(Sigma0) != ncol(Sigma0) ) {
        stop('Sigma0 must be square')
      }
      if ( any(eigen(Sigma0)$values <= 0) ) {
        stop('Sigma0 must be positive definite')
      }
      if ( length(formula.list) != nrow(Sigma0) ) {
        stop('Sigma0 and formula.list must have the same dimensions')
      }
    }
    
    ## Check 0 < a0 <= 1
    if ( a0 <= 0 | a0 > 1 ) { stop('a0 must be larger than 0 and less than or equal to 1') }
  
    ## check burnin is non-negative integer
    if ( !( burnin >= 0 & burnin %% 1 == 0 ) ) {
      stop('burnin must be a non-negative integer')
    }
  
    if ( !(thin >= 1 & burnin %% 1 == 0 ) ) {
      stop('burnin must be a positive integer')
    }
  
  
  ## Number of endpoints / obs
  n <- nrow(data)
  J <- length(formula.list)
  
  ## Helper functions to get Xlist and Y matrix
  get_Xlist <- function(formula.list, data) {
    Xlist = lapply(formula.list, function(f, data) { model.matrix(f, data = data) }, data = data )
  }
  get_Ymat <- function(formula.list, data) {
    ynames <- sapply(formula.list, function(f, data) all.vars(f)[1])
    Ymat <- as.matrix(data[, ynames])
  }
  
  ## Obtain list of X's and extract names
  Xlist0 <- get_Xlist(formula.list, histdata)
  Xlist <- get_Xlist(formula.list, data)
  ynames <- sapply(formula.list, function(f) all.vars(f)[1] )       ## names of responses
  

  ## Get names of predictors as j.xname as vector
  pvec <- sapply(Xlist, function(x) ncol(x))
  J <- length(formula.list)
  Xnames <- lapply(formula.list, function(f, data) colnames(model.matrix(f, data)), data = data)
  for(i in 1:J) {
    Xnames[[i]] <- paste0(i, '.', Xnames[[i]])
  }
  Xnames <- unlist(Xnames)
  
  
  ## Get response matrix (Y1, ..., YJ)
  Y0 <- get_Ymat(formula.list, histdata)
  Y <- get_Ymat(formula.list, data)
  y0 <- as.vector(Y0)
  y <- as.vector(Y)
  
  X0 <- as.matrix(Matrix::bdiag(Xlist0))
  X <- as.matrix(Matrix::bdiag(Xlist))
  XtX <- crossprod(rlist::list.cbind(Xlist))
  X0tX0 <- crossprod(rlist::list.cbind(Xlist0))
  
  if ( is.null(Sigma0) ) {
    resmat <- sapply(formula.list, function(f, data) stats::resid(lm(f, data) ), data = data )
    Sigma0 <- crossprod(resmat) / nrow(data)
  }
  
  sample <- sur_sample_gibbs_cpp( Sigma0, M, X, X0, XtX, X0tX0, Y, Y0, y, y0, a0, pvec, burnin, thin )
  
  betasample <- sample$betasample
  colnames(betasample) <- Xnames
  Sigmalist <- sample$Sigmalist
  
  
  res <- list(
    'betadraw' = betasample,
    'Sigmalist' = Sigmalist,
    'formula.list' = formula.list,
    'ynames' = ynames,
    'Xnames' = Xnames,
    'num.covars' = pvec,
    'n' = n,
    'J' = J,
    'method' = 'gibbs'
  )
  
  ## Set the name for the class
  class(res) <- append(class(res), "surbayes" )
  return(res)
}

