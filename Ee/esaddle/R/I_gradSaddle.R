#####
# Computes the gradient of the saddlepoint density
# @description Computes the gradient of the saddlepoint density
#
# @param lambda Point at which the CGF is evaluated (d-dimensional vector).
# @param X (n by d) matrix containing the data.
# @param deriv If TRUE the gradient of the empitical CGF wrt y (and at y) is returned.
#              Otherwise the values of the empirical CGF (and possibly of its derivatives wrt
#              lambda) at lambda is returned.
# @param onlyDlamDy if TRUE only dLambda/dY is computed.
# @param addList = list of additional (optional) arguments: 
#         \itemize{
#         \item{ \code{invCOV} }{The inverse of kum2;}
#         \item{ \code{y} }{The point at which the underlying empirical saddlepoint is evaluated;}
#         \item{ \code{grad} }{The decay rate of the saddlepoint. See ?dsaddle for details;}
#         }
# @return The gradient of the saddlepoint. 
# @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#
.gradSaddle <- cmpfun( function(y, lambda, X, decay, extra, onlyDlamDy = FALSE, mixMethod = "mse") {
  
  if(!is.vector(lambda)) lambda <- as.vector(lambda) 
  if (!is.matrix(X)) X <- matrix(X, length(X), 1)
  n <- nrow(X)
  d <- ncol(X)
  
  stopifnot(d == length(lambda))
  
  kum1 <- rep(0, d)
  kum2 <- diag(1, d)
  
  # Mixture and derivative of mixture function wrt y
  tmp <- .ecgfMix(y, decay = decay, method = mixMethod, deriv = TRUE, m = n)
  mix <- tmp[["mix"]]
  DmixDy <- tmp[["DmixDy"]]
  
  elx <- drop(extra$elx) # exp(lambda'x_i - alpha) vector
  
  # Original K, K' and K'', before tilting with Gaussian cgf
  tmp_K <- extra$tmp_K
  tmp_dK <- drop( extra$tmp_dK )
  tmp_d2K <- extra$tmp_d2K
  
  # Scaling K''
  dd <- extra$dd
  D <- diag(dd, nrow = d, ncol = d)
  DD <- extra$DD
  d2KQR <- extra$d2KQR
  
  # Derivative of lambda wrt y
  DlamDy <- dd * qr.solve(d2KQR, dd * (diag(1, d) - tcrossprod(tmp_dK - lambda, DmixDy)), tol = 0)
  
  if( onlyDlamDy ) return( DlamDy )
  
  # Computing K'''
  d3K <- array(NA, c(d, d, d) )
  
  for(ff in 1:d)
  {
    A <- crossprod(X, (elx*X)*X[ , ff] ) / sum(elx)
    B <- -( tmp_d2K + tcrossprod(tmp_dK, tmp_dK) ) * tmp_dK[ff] 
    C <- - tcrossprod(tmp_d2K[ , ff], tmp_dK) - t( tcrossprod(tmp_d2K[ , ff], tmp_dK) )
    d3K[ , , ff] <- A + B + C 
  }  
  
  # Computing derivatives of K'' wrt y and the gradient of the saddlepoint log-density wrt y
  d2Kdy <- matrix(0, d, d)
  DsadDy <- numeric(d)
  for( ff in 1:d )
  {
    d2Kdy <- d2Kdy * 0 
    for(zz in 1:d)
    {
      d2Kdy <- d2Kdy + d3K[ , , zz] * DlamDy[zz, ff] # Should it be [ff, zz] ??
    }
    
    DsadDy[ff] <- - mix * 0.5 * .Trace( D %*% qr.solve(d2KQR, D %*% d2Kdy, tol = 0) ) 
  }
  
  DsadDy <- DsadDy - 
    lambda + 
    DlamDy %*% (extra$dK - y) -  
    0.5 * .Trace( D %*% qr.solve(d2KQR, D %*% ( tmp_d2K - kum2 ), tol = 0) ) * DmixDy + 
    DmixDy * drop( tmp_K - crossprod(kum1, lambda) - 0.5 * crossprod(lambda, kum2%*%lambda) )
  
  return( list("DsadDy" = DsadDy, "DLamDy" = DlamDy, "DmixDy" = DmixDy) )
  
})
