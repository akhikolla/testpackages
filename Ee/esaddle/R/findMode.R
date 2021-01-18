#####
#' Finding the mode of the empirical saddlepoint density
#' @description Given a sample from a d-dimensional distribution, the routine
#'              finds the mode of the corresponding Extended Empirical Saddlepoint (EES) density.
#' @param X an n by d matrix containing the data.
#' @param init d-dimensional vector containing the starting point for the optimization. By default
#'             it is equal to \code{colMeans(X)}.
#' @param decay rate at which the SPA falls back on a normal density. Should be a positive number. See Fasiolo et al. (2016)
#'              for details.
#' @param method optimization method used by \code{stats::optim()}, see ?optim for details. By default it is "BFGS".
#' @param hess if TRUE also an estimate of the Hessian at the mode will be returned.
#' @param sadControl list corresponding to the \code{control} argument in the \code{dsaddle} function.
#' @param ... Extra arguments to be passed to the optimization routine \code{stats::optim}. 
#' @return A list where \code{mode} is the location of mode of the empirical saddlepoint density,
#'         \code{logDens} is the log-density at the mode and \code{hess} (present only if argument \code{hess==TRUE})
#'         is the approximate Hessian at convergence. The other entries are the same as for \code{stats::optim}.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @references Fasiolo, M., Wood, S. N., Hartig, F. and Bravington, M. V. (2016). 
#'             An Extended Empirical Saddlepoint Approximation for Intractable Likelihoods. ArXiv http://arxiv.org/abs/1601.01849.
#' @examples 
#' # library(esaddle)
#' set.seed(4141)
#' x <- rgamma(1000, 2, 1)
#' 
#' # Fixing tuning parameter of EES
#' decay <-  0.05
#' 
#' # Evaluating EES at several point
#' xSeq <- seq(-2, 8, length.out = 200)
#' tmp <- dsaddle(y = xSeq, X = x, decay = decay, log = TRUE)  # Un-normalized EES
#' 
#' # Plotting true density, EES and normal approximation
#' plot(xSeq, exp(tmp$llk), type = 'l', ylab = "Density", xlab = "x")
#' lines(xSeq, dgamma(xSeq, 2, 1), col = 3)
#' suppressWarnings( rug(x) )
#' legend("topright", c("EES", "Truth"), col = c(1, 3), lty = 1)
#' 
#' # Find mode and plot it
#' res <- findMode(x, init = mean(x), decay = decay)$mode
#' abline(v = res, lty = 2, lwd = 1.5)
#' @export
#'
findMode <- function(X, decay, init = NULL, method = "BFGS", hess = FALSE, sadControl = list(), ...)
{
  switch(class(X)[1],
         "matrix"  = theData <- X,
         "numeric" = theData <- matrix(X, length(X), 1),
         stop("X should be either of class \"vector\" or a matrix")
  )
  
  nDims <- ncol(theData)
  nObs <- nrow(theData)
  
  if( nDims > nObs )
  {
    warning("The input data should be seems to have more dimensions (columns) than data points (rows): I'm transposing it!")
    theData <- t(theData)
  }
  
  if(is.null(init)){ 
    init <- colMeans(theData) 
  } else { 
    stopifnot( is.vector(init), length(init) == ncol(theData) ) 
  }
  
  objFun  <- function(x) -dsaddle(y = as.numeric(x), X = theData, control = sadControl, decay = decay, log = TRUE)$llk
  
  objGrad <- function(x) -dsaddle(y = as.numeric(x), X = theData, control = sadControl, decay = decay, deriv = TRUE)$grad
  
  optOut <- optim(par = init, fn = objFun, gr = objGrad, method = method, ...)
      
  names(optOut)[c(1, 2, 3)] <- c("mode", "logDens", "numEval")
  
  optOut$logDens <- - optOut$logDens
  if(hess)
  {
    optOut$hess <- - .hessFromGrad(optOut$mode, objGrad, eps = sqrt(.Machine$double.eps), ...)
  }
  
  return(optOut)
}