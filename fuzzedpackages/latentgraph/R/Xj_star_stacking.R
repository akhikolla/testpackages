Xj_star_stacking <- function(dat,n,R,p){ #stack the data by replicates
  for(i in 1:n){
    dat[R*i,] <- 0 #put 0 for the 1st replicate of every R replicates
  }
  a <- matrix(0,1,p) #1st row is 0
  b <- rbind(a,dat)
  Xj_star <- b[-dim(b)[1],] #exclude the last row
  return(Xj_star)
}

