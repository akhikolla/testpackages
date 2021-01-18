#'@method predict PUfit
#'@export
predict.PUfit <- function(object,newdata,lambda=NULL,type=c("response","link"),...){
  
  type <- match.arg(type)
  coef <- coef(object,lambda=lambda)
  eta <- cbind(rep(1,nrow(newdata)),newdata)%*%coef
  p <- matrix(NA,ncol=ncol(eta),nrow=nrow(eta))
  for (i in 1:ncol(eta)){
    p[,i] = 1/(1+exp(-eta[,i]))
  }
  
  if(type=="link"){
    return(eta)
  } else if(type=="response"){
    return(p)
  } else{
    stop("wrong type")
  } 
}
#'@export
predict.cvPUfit <- function(object,newdata,lambda=NULL,type=c("response","link"),...){
  
  type <- match.arg(type)
  coef <- coef(object,lambda=lambda)
  eta <- cbind(rep(1,nrow(newdata)),newdata)%*%coef
  p <- matrix(NA,ncol=ncol(eta),nrow=nrow(eta))
  for (i in 1:ncol(eta)){
    p[,i] = 1/(1+exp(-eta[,i]))
  }
  
  if(type=="link"){
    return(eta)
  } else if(type=="response"){
    return(p)
  } else{
    stop("wrong type")
  } 
}
