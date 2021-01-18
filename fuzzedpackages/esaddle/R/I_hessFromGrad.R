#### Calculates Hessian, given a function that calculates the gradient
# INPUT
# x: point at which the gradient is calculated
# gr: gradient function
# eps: step size

.hessFromGrad <- function(x, gr, eps = sqrt(.Machine$double.eps), ...)
{
  d <- length(x)
  
  H <- matrix(NA, d, d)
  p <- numeric(d)
  
  for(ii in 1:d)
  {
    p[ii] <- eps
    
    H[ , ii] <- ( (gr(x + p, ...) - gr(x - p, ...)) / (2*eps) )
    
    p[ii] <- 0
  }
  
  H <- (H + t(H)) / 2
  
  return( H )
  
}
