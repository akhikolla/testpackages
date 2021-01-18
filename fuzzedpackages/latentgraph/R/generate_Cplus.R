generate_Cplus <- function(R,n){
  C <- generate_C(R,n)
  A <- matrix(1, nrow=1,ncol=R)
  E <- matrix(0, nrow=n,ncol=n*R)
  for(i in 1:n){
    E[i, (1+(i-1)*R):(R+(i-1)*R)] = A #put 1's for every T columns, row by row
  }

  C_tilde <- rbind(C, E)
  C_plus <- pinv(C_tilde) #computes the pseudoinverse

  return(C_plus)
}
