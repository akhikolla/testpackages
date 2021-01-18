updateS <- function(R,L,Gamma,lambda1,rho){
  A <- R+L+Gamma/rho
  B <- lambda1/rho
  S <- Soft(A,B)
  return(S)
}
