# Calculates the laplace approximation given the value (f) of a function at its
# maximum and the Hessian (H)

.laplApprox <- function(f, H, log = FALSE)
{
  
  d <- ncol(H)
  
  out <- d/2 * log(2*pi) - determinant(H, logarithm = T)$modulus/2 + f
  
  if( !log )
  {
    out <- exp( out )
  }
  
  return( as.numeric(out) )
  
}