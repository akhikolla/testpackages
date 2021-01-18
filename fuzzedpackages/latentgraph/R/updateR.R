updateR <- function(S,L,Sigmahat,Gamma,rho){
  C <-  S-L-(Sigmahat+Gamma)/rho
  a <- eigen(C)
  D <- diag(a$values)
  U <- a$vectors
  R <- 1/2*U%*%(D+sqrt(D*D+4/rho*diag(rep(1,nrow(Sigmahat)))))%*%t(U)
  return(R)
}
