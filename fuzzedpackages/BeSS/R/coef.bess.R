coef.bess=function(object, sparse=TRUE, type = c("ALL", "AIC", "BIC", "EBIC"),...)
{
  type <- match.arg(type)
  if(!is.null(object$coef0)){
    beta=rbind(intercept=object$coef0, object$beta)
    rownames(beta)[1] = "(intercept)"
  } else beta=object$beta
   if(sparse==TRUE)
   {
     beta=Matrix(beta,sparse = TRUE)
     if(type == "ALL") {
     return(beta)
     }else return(Matrix(beta[,which.min(object[[type]])], sparse = TRUE))
   }else {
     if(type == "ALL") {
       return(beta)
     }else return(beta[,which.min(object[[type]])])
   }
}

coef.bess.one=function(object,sparse=TRUE, ...)
{
  if(!is.null(object$coef0)){
    beta=c(intercept=object$coef0, object$beta)
    names(beta)[1] = "(intercept)"
    } else beta=object$beta
  if(sparse==TRUE)
  {
    beta=matrix(beta,byrow =TRUE)
    beta=Matrix(beta,sparse = TRUE)
    return(beta)
  }else return(beta)
}
