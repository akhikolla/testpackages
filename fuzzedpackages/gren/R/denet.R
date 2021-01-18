# density function of the elastic net prior distribution
denet <- function(x, lambda1=1, lambda2=1, log=FALSE) {
  if(log) {
    dens <- 0.5*log(lambda2) - log(2) - 0.5*lambda1*abs(x) - 0.5*lambda2*x^2 +
      dnorm(0.5*lambda1/sqrt(lambda2), log=TRUE) - 
      pnorm(-0.5*lambda1/sqrt(lambda2), log.p=TRUE)
  } else {
    dens <- 0.5*sqrt(lambda2)*exp(-0.5*(lambda1*abs(x) + lambda2*x^2))*
      exp(dnorm(0.5*lambda1/sqrt(lambda2), log=TRUE) - 
            pnorm(-0.5*lambda1/sqrt(lambda2), log.p=TRUE))
  }
  return(dens)
}