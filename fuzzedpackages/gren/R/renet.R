# sampling function of the elastic net prior distribution
renet <- function(n, lambda1=1, lambda2=1) {
  if(lambda1==0) {
    beta <- rnorm(n, 0, sqrt(1/lambda2))
  } else if(lambda2==0) {
    u <- runif(n) - 0.5
    beta <- -2*sign(u)*log(1 - 2*abs(u))/lambda1
  } else {
    tau <- rtau(n, lambda1, lambda2)
    sigma <- (tau - 1)/(lambda2*tau)
    beta <- rnorm(n, 0, sqrt(sigma))
  }
  return(beta)
}