semigraph_GGM <- function(newdata,n,p,lambda=NULL){

  if(length(lambda)==1){
    lambda <- rep(lambda,p)
  }

  beta <- NULL
  if(is.null(lambda)==TRUE){
    lambda <- rep(0,p)
  }

  for(j in 1:p){

    restemp <- glmnet(newdata[[j]][[2]],newdata[[j]][[1]],family=c("binomial"),lambda=lambda[j],standardize=FALSE,intercept=FALSE)

    betatemp <- coef(restemp)[2:(p)]
    if(j==1){
      betatemp <- c(0,betatemp)
    }
    if(j!=1 && j < p){
      betatemp <- c(betatemp[1:(j-1)],0,betatemp[(j):(p-1)])
    }
    if(j==p){
      betatemp <- c(betatemp,0)
    }
    beta <- cbind(beta,betatemp)

  }
  colnames(beta) <- paste("feature",c(1:p),sep=" ")

  return(list(beta=beta,lambda=lambda))

}
