generate_C <- function(R,n){ #calculate the distance between different thetas...minimize to make piecewise representation
  D <- matrix(0,nrow=n*R,ncol=n*R)

  for(i in 1:(n*R-1)){
    D[i,i] <- -1
    D[i,i+1] <- 1
  }

  for(i in 1:n){
    D[i*R,] <- NA
  }
  C <- na.omit(D) #omit the lines that are the 1st of every T replicates

  return(C)
}
