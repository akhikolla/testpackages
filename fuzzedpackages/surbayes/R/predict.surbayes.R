#' Get predictive posterior samples
#'
#' This function returns a list of new data sets by sampling
#' from the posterior predictive density of Y | Y0, Xnew.
#' 
#' @importFrom stats as.formula lm model.matrix
#' @param object Result from calling \code{sur_sample}
#' @param newdata \code{data.frame} of new data
#' @param nsims number of posterior draws to take. The default and minimum is 1. The maximum is the number of simulations in surbayes
#' @param ... further arguments passed to or from other methods
#' @return \code{n x J x nsims} \code{array} of predicted values
#' 
#' @examples 
#' 
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
#' ## Fit model
#' out = sur_sample( formula.list, data, M = M )
#' 
#' ## Obtain predictions
#' pred = predict(out, data, nsims = 1)
#' 
#' @export
predict.surbayes <- function( object, newdata, nsims = 1, ... ) {
  
  if(!'surbayes' %in% class(object)) {
    stop("Must provide object of type surbayes")
  }
  
  if ( nsims < 1 | nsims > nrow(object$betadraw ) ) {
    stop("Number of posterior samples must be between 1 and number of simulations in object")
  }
  
  if(class(newdata) != 'data.frame') {
    stop("newdata must be a data.frame")
  }
  
  ## Obtain new design matrices
  Xnewlist <- lapply( object$formula.list, model.matrix, data = newdata )
  Xnew <- Matrix::bdiag(Xnewlist)
  
  ## Obtain matrix of means; each column is one draw
  M <- as.matrix( Xnew %*% t( object$betadraw[1:nsims, , drop = FALSE] ) ) 
  
  ## Get various statistics to pass onto C++ function
  n <- nrow(Xnew)
  J <- object$J
  
  ## Call C++ helper function to loop through samples and draw from historical posterior
  res <- predict_surbayes_cpp( M, object$Sigmalist, n, J, nsims )
  dimnames(res)[[2]] <- object$ynames
  
  return(res)
}
