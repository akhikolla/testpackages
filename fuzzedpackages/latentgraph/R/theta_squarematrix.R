theta_squarematrix <- function(y,p,rule = ""){ #y is the theta matrix; d is an empty(0's) square matrix
  d <- matrix(0, nrow = p, ncol = p)
  d[,1] <- c(0,y[,1])
  for(j in 2:p){
    d[,j] <- append(y[,j],0,after = j-1) #get 0's on the diagonal b/c we don't care about relationship with itself
  }

  return(d)
}
