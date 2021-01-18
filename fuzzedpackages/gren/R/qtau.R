# quantile function of the truncated gamma with shape=1/2, 
# scale=8*lambda2/lambda1^2, and domain=(1,infinity)
qtau <- function(p, lambda1, lambda2, log.p=FALSE) {
  scale <- 8*lambda2/lambda1^2
  if(log.p) {
    aux <- pnorm(-sqrt(2/scale), log.p=log.p)
    tau <- 0.5*scale*qnorm(exp(aux) - exp(aux + p))^2
  } else {
    tau <- 0.5*scale*qnorm((1 - p)*pnorm(-sqrt(2/scale)))^2
  }
  return(tau)
}