transformEachIndex <- function(X,j,n,p,R){
  if(length(R)<=1)
    stop("Error! R is an n-vector indicating the number of replicates for each subject")
  finaly <- NULL
  finalZ <- NULL

  for(i in 1:n){
    ithobs <- X[[i]]

    # Get the pairwise difference
    temp <- transformEachObs(j,p,R[i],ithobs)
    finaly <- c(finaly,temp$ytilde)
    finalZ <- rbind(finalZ,temp$Xtilde)
  }

  return(list(y=finaly,Z=finalZ))
}
