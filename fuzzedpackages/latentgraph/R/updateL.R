updateL <- function(R,S,Gamma,lambda2,rho){
  C <- S-R-Gamma/rho
  a <- eigen(C)
  D <- diag(a$values)
  U <- a$vectors
  pmax(D-lambda2/rho,0)
  L <- U%*%(pmax(D-lambda2/rho,0))%*%t(U)
  return(L)
}
