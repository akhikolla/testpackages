# sampling function of the truncated gamma with shape=1/2, 
# scale=8*lambda2/lambda1^2, and domain=(1,infinity)
rtau <- function(n, lambda1, lambda2) {
  u <- runif(n)
  tau <- qtau(u, lambda1, lambda2)
  return(tau)
}