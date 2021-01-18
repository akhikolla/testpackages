transformDataReplicates <- function(X,n,p,R){
  if(length(R)<=1)
    stop("Error! R is an n-vector indicating the number of replicates for each subject")
  newdata <- vector('list',p)
  for(j in 1:p){
    chdata <- transformEachIndex(X,j,n,p,R)
    newdata[[j]] <- vector('list',2)
    newdata[[j]][[1]] <- chdata$y
    newdata[[j]][[2]] <- chdata$Z
  }
  return(newdata)
}
