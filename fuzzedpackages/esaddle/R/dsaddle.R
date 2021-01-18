#####
#' Evaluating the Extended Empirical Saddlepoint (EES) density 
#' @description Gives a pointwise evaluation of the EES density (and optionally of its gradient) at one or more 
#'              locations.
#'
#' @param y points at which the EES is evaluated (d dimensional vector) or an n by d matrix, each row indicating
#'          a different position.
#' @param X n by d matrix containing the data.
#' @param decay rate at which the EES falls back on a normal density approximation, fitted to \code{X}. 
#'              It must be a positive number, and it is inversely proportional to the complexity of the fit.
#'              Setting it to \code{Inf} leads to a Gaussian fit.  
#' @param deriv If TRUE also the gradient of the log-saddlepoint density is returned.
#' @param log If TRUE the log of the saddlepoint density is returned.
#' @param normalize If TRUE the normalizing constant of the EES density will be computed. FALSE by 
#'                  default.
#' @param control A list of control parameters with entries:
#'         \itemize{
#'         \item{ \code{method} }{the method used to calculate the normalizing constant. 
#'                                Either "LAP" (laplace approximation) or "IS" (importance sampling).}
#'         \item{ \code{nNorm} }{if control$method == "IS", this is the number of importance samples used.}
#'         \item{ \code{tol} }{the tolerance used to assess the convergence of the solution to the saddlepoint equation.
#'                             The default is 1e-6.}
#'         \item{ \code{maxit} }{maximal number of iterations used to solve the saddlepoint equation.
#'                               The default is 100;}
#'         \item{ \code{ml} }{Relevant only if \code{control$method=="IS"}. n random variables are generated from 
#'                            a Gaussian importance density with covariance matrix \code{ml*cov(X)}. 
#'                            By default the inflation factor is \code{ml=2}.}
#'         }
#' @param multicore  if TRUE the empirical saddlepoint density at each row of y will be evaluated in parallel.
#' @param ncores   number of cores to be used.
#' @param cluster an object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster. 
#' @return A list with entries:
#'         \itemize{
#'         \item{ \code{llk} }{the value of the EES log-density at each location y;}
#'         \item{ \code{mix} }{for each location y, the fraction of saddlepoint used: 
#'                             1 means that only ESS is used and 0 means that only a Gaussian fit is used;}
#'         \item{ \code{iter} }{for each location y, the number of iteration needed to solve the 
#'                              saddlepoint equation;}
#'         \item{ \code{lambda} }{an n by d matrix, where the i-th row is the solution of the saddlepoint 
#'                                equation corresponding to the i-th row of y;}
#'         \item{ \code{grad} }{the gradient of the log-density at y (optional);}
#'         \item{ \code{logNorm} }{the estimated log normalizing constant (optional);}
#'         }
#' @references Fasiolo, M., Wood, S. N., Hartig, F. and Bravington, M. V. (2016). 
#'             An Extended Empirical Saddlepoint Approximation for Intractable Likelihoods. ArXiv http://arxiv.org/abs/1601.01849.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon N. Wood.
#' @examples 
#' library(esaddle)
#' 
#' ### Simple univariate example
#' set.seed(4141)
#' x <- rgamma(1000, 2, 1)
#' 
#' # Evaluating EES at several point
#' xSeq <- seq(-2, 8, length.out = 200)
#' tmp <- dsaddle(y = xSeq, X = x, decay = 0.05, log = TRUE)  # Un-normalized EES
#' tmp2 <- dsaddle(y = xSeq, X = x, decay = 0.05,             # EES normalized by importance sampling
#'                 normalize = TRUE, control = list("method" = "IS", nNorm = 500), log = TRUE)
#' 
#' # Plotting true density, EES and normal approximation
#' plot(xSeq, exp(tmp$llk), type = 'l', ylab = "Density", xlab = "x")
#' lines(xSeq, dgamma(xSeq, 2, 1), col = 3)
#' lines(xSeq, dnorm(xSeq, mean(x), sd(x)), col = 2)
#' lines(xSeq, exp(tmp2$llk), col = 4)
#' suppressWarnings( rug(x) )
#' legend("topright", c("EES un-norm", "EES normalized", "Truth", "Gaussian"), 
#'         col = c(1, 4, 3, 2), lty = 1)
#' @export
#'
dsaddle <- function(y, X,  decay, deriv = FALSE, log = FALSE, 
                    normalize = FALSE, control = list(), 
                    multicore = !is.null(cluster), ncores = detectCores() - 1, cluster = NULL) {
  ## X[i,j] is ith rep of jth variable; y is vector of variables.
  ## evaluate saddle point approximation based on empirical CGF
  if( !is.matrix(X) ) X <- matrix(X, length(X), 1)
  
  d <- ncol(X)
  
  if( !is.matrix(y) ){ 
    if(d == 1){
      y <- as.matrix(y)
    } else {
      if(length(y) == d) y <- matrix(y, 1, d) else stop("y should be a matrix n by d.")
    }
  }
  
  ny <- nrow( y )
  
  # Offsetting dimensionality, so decay stays pretty much at the same level for any d.
  decayI <- decay / ( d ^ 2 )
  
  # Setting up control parameter
  ctrl <- list( "method" = "LAP", 
                "nNorm" = 100 * ncol(X), 
                "tol" = 1e-6,
                "maxit" = 100,
                "ml" = 2)
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control, verbose = FALSE)
  
  if( multicore ){ 
    # Force evaluation of everything in the environment, so it will available on cluster
    .forceEval(ALL = TRUE)
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "esaddle", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoParallel(cluster)
  }
  
  # Pre-calculating covariance and normalizing
  #iCov <- .robCov(t(X), alpha2 = 4, beta2 = 1.25)
  iCov <- .robCov(t(X), alpha = 10, alpha2 = 10, beta2 = 1.25)
  
  # Saving originals
  iX <- X
  iy <- y
  
  # Weighting the statistics in order to downweight outliers
  X <- iCov$weights * X
  
  # Creating normalized version
  y <- t( iCov$E %*% (t(y) - iCov$mY) )
  X <- t( iCov$E %*% (t(X) - iCov$mY) )
  
  # If there are some statistics with zero variance we remove them
  if( length(iCov$lowVar) )
  {
    stop("The columns of X indexed", iCov$lowVar, "have zero variance.")
    #y <- y[-iCov$lowVar]
    #X <- X[ , -iCov$lowVar, drop = FALSE]
  }
  
  out <- list()
  # Divide saddlepoint evaluations between cores
  withCallingHandlers({
    out <- alply(y, 1, .dsaddle, .parallel = multicore,
                 # Args for .dsaddle()
                 X = X,
                 tol = ctrl$tol,
                 maxit = ctrl$maxit,
                 decay = decayI,
                 deriv = deriv)}, warning = function(w) {
                   # There is a bug in plyr concerning a useless warning about "..."
                   if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
                     invokeRestart("muffleWarning")
                 })
  
  # Close the cluster if it was opened inside this function
  if(multicore && clusterCreated) stopCluster(cluster)
  
  
  
  # Create output list
  out <- list( "llk" = sapply(out, "[[", "llk"),
               "mix" = sapply(out, "[[", "mix"),
               "niter" = sapply(out, "[[", "niter"),
               "lambda" = t(sapply(out, "[[", "lambda")),
               "grad" = if(deriv) t( sapply(out, "[[", "DsadDy") ) else NULL)
  
  # Correct the log-likelihood for the missing normalizing constant
  if( normalize ){
    
    stopifnot( (ctrl$method %in% c("IS", "LAP")) )
    
    # Log-normalizing constant by importance sampling
    if(ctrl$method == "IS") 
    {
      aux <- rmvn(ctrl$nNorm, iCov$mY, ctrl$ml*iCov$COV)
      
      logNorm <- log( .meanExpTrick( 
        dsaddle(y = aux, X = iX, decay = decay, deriv = FALSE, 
                log = TRUE, normalize = FALSE, control = ctrl, 
                multicore = multicore, ncores = ncores, cluster = cluster)$llk - dmvn(aux, iCov$mY, ctrl$ml*iCov$COV, log = TRUE) )
      )
      
    }
    
    # Log-normalizing constant by Laplace approximation
    if(ctrl$method == "LAP")
    {
      tmp <- findMode(X = iX, init = iCov$mY, decay = decay, 
                      sadControl = list("tol" = ctrl$tol, "maxit" = ctrl$maxit), hess = T)
      logNorm <- .laplApprox(tmp$logDens, tmp$hess, log = TRUE)
    }
    
    out$logNorm <- logNorm
    out$llk <- out$llk - logNorm
    
  }
  
  # Adjusting log-lik and its gradient to correct for normalization
  out$llk <- out$llk + sum(log(abs(diag(iCov$E))))
  if( deriv ) out$grad <- drop( drop(out$grad) %*% iCov$E )
  
  if( !log ) out$llk <- exp( out$llk )
  
  return( out )
  
}




##########
# INTERNAL
##########

.dsaddle <- cmpfun(function(y, X, tol, decay, 
                            deriv = FALSE, mixMethod = "mse", 
                            maxit = 100, lambda = NULL) {
  ## X[i,j] is ith rep of jth variable; y is vector of variables.
  ## evaluate saddle point approximation based on empirical CGF
  
  if( !is.vector(y) ) y <- as.vector(y)
  
  d <- length(y)
  
  if( !is.matrix(X) ){
    if(d > 1){ 
      stop("Error: simulated data must be entered in matrix form")
    }else{
      X <- matrix(X, length(X), 1)
    }
  }
  
  n <- nrow(X)
  
  # Initial guess of the root is the solution to the Gaussian case
  # the gain is one step less of Newton on average.
  if( is.null(lambda) ) lambda <- y
  
  mu <- rep(0, d)
  sig <- diag(1, d)
  
  # Choose the mixture of saddlepoint-normal, mix \in [0, 1]
  mix <- .ecgfMix(y, decay = decay, method = mixMethod, deriv = FALSE, m = d)$mix
  
  b <- .ecgf(lambda, X, kum1 = mu, kum2 = sig, mix = mix, grad = 2)
  
  ## Newton loop to minimize K(lambda) - t(lambda)%*%y or solve dK(lambda) = y wrt lambda
  kk <- jj <- 0
  # Convergence test: see [con_test] below.
  while( any( abs(b$dK-y) > tol ) && kk < maxit ) 
  { 
    kk <- kk + 1
    
    # Build scaling vector and matrix
    dd <- diag(b$d2K)^-0.5
    DD <- tcrossprod(dd, dd)
    
    d2KQR <- qr(DD * b$d2K, tol = 0)
    
    # Try solve scaled linear system fast, if that doesn't work use QR decomposition.
    d.lambda <- -  dd * drop(qr.solve(d2KQR, dd*(b$dK-y), tol = 0))
    
    lambda1 <- lambda + d.lambda ## trial lambda
    
    b1 <- .ecgf(lambda1, X, kum1 = mu, kum2 = sig, mix = mix, grad = 2)
    if ( sum( abs(b1$d2K) ) == 0 ) return(NA) ## spa breakdown (possibly too tight)
    
    jj <- 1
    c1 <- 10^-4
    alpha <- 1
    rho <- 0.5
    ## Line search checking Arminjo condition and step halving
    while( ( b1$K - crossprod(lambda1, y) ) > 
           ( b$K - crossprod(lambda, y) ) + c1 * alpha * crossprod(d.lambda, drop(b$dK-y)) && jj < 50)  
    {
      jj <- jj + 1
      alpha <- alpha * rho
      d.lambda <- d.lambda * alpha
      lambda1 <- lambda + d.lambda
      b1 <- .ecgf(lambda1, X, kum1 = mu, kum2 = sig, mix = mix, grad = 2)
    }
    
    ## now update lambda, K, dK, d2K
    lambda <- lambda1 
    b <- b1
    
  } ## end of Newton loop
  ## lambda is the EES lambda...
  
  if(kk > 50 || jj > 20) warning(paste("The convergence of the saddlepoint root-finding is quite slow! \n",
                                       "Outer root-finding Newton-Raphson:", kk, "iter \n",
                                       "Inner line search for Arminjo condition:", jj, "iter"))
  
  # We need to recompute the QR decomposition
  # Build scaling vector and matrix
  dd <- diag(b$d2K)^-0.5
  DD <- tcrossprod(dd, dd)
  d2KQR <- qr(DD * b$d2K, tol = 0)
  
  # Determinant of b$d2K from its qr decomposition
  logDet <- sum(log(abs(diag(qr.R(d2KQR))))) - 2 * sum(log(dd))
  
  spa <- b$K - crossprod(lambda, y) - 0.5 * log( 2*pi ) * d - 0.5 * logDet
  
  # Additional stuff needed by .gradSaddle()
  b[ c("dd", "DD", "d2KQR") ] <- list(dd, DD, d2KQR)
  
  out <- list("llk" = drop(spa), "mix" = mix, "niter" = kk, "lambda" = lambda, "extra" = b)
  
  if(deriv)
  {
    
    tmp <- .gradSaddle(y = y, lambda = lambda, X = X, decay = decay, mixMethod = mixMethod, extra = b)
    
    out <- c(out, tmp)
    
  }
  
  return( out )
  
})