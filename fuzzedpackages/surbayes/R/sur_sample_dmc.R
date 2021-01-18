#' Sample SUR model via direct Monte Carlo
#'
#' This function samples from the posterior of a SUR model using the DMC method of Ando and Zellner (2010)
#' 
#' @importFrom stats as.formula lm model.matrix terms
#' @param formula.list A list of formulas, each element giving the formula for the corresponding endpoint.
#' @param data A \code{data.frame} containing all the variables contained in \code{formula.list]}
#' @param M Number of samples to be drawn
#' @return A list. First element is posterior draws. Second element is list of JxJ covariance matrices. Other elements are helpful statistics about the SUR model to pass to other functions.
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
#' ## fit using historical data as current data set--never done in practice
#' out = sur_sample_powerprior( formula.list, data, histdata = data, M = M )
#' 
#' @import Rcpp
#' @importFrom Matrix bdiag
#' @importFrom rlist list.cbind
#' @export
sur_sample_dmc <- function(
  formula.list,
  data,
  M = 1000
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
    
    ## Check if M is positive integer
    if ( M <= 0 | M %% 1 != 0 ) { stop('M must be a positive integer') }
  
  
  
  ## Get various needed parameters to pass to functions
  Xlist <- lapply(formula.list, model.matrix, data = data)
  Xnames.list <- lapply(formula.list, function(f) all.vars(f)[-1])  ## names of predictors
  ynames <- sapply(formula.list, function(f) all.vars(f)[1] )       ## names of responses
  Y <- as.matrix(data[, ynames])                    ## matrix of responses (Y1, ..., YJ)
  y <- as.vector(Y)                                 ## vector (y1', ..., yJ')'
  X <- as.matrix(Matrix::bdiag(Xlist))              ## full design matrix
  XtX <- crossprod(rlist::list.cbind(Xlist))        ## p x p matrix of crossproducts
  p.j <- sapply(Xlist, ncol)                        ## number of covars / endpoint
  n <- nrow(data)                                   ## number of observations
  J <- ncol(Y)                                      ## number of endpoints
  
  ## Get names of predictors as j.xname as vector
  pvec <- sapply(Xlist, function(x) ncol(x))
  J <- length(formula.list)
  Xnames <- lapply(formula.list, function(f, data) colnames(model.matrix(f, data)), data = data)
  for(i in 1:J) {
    Xnames[[i]] <- paste0(i, '.', Xnames[[i]])
  }
  Xnames <- unlist(Xnames)
  
  
  ## Call C++ function to sample beta, Sigma | y, X
  sample.list <- sur_sample_cpp(Y, Xlist, y, X, XtX, p.j, M)
  
  ## Name the columns of beta samples
  colnames(sample.list$betadraw) <- Xnames
  
  res <- list(
    'betadraw' = sample.list$betadraw,
    'Sigmalist' = sample.list$covlist,
    'formula.list' = formula.list,
    'ynames' = ynames,
    'Xnames' = Xnames,
    'num.covars' = p.j,
    'n' = n,
    'J' = J,
    'method' = 'dmc'
  )
  ## Set the name for the class
  class(res) <- append(class(res),"surbayes")
  return(res)
}
