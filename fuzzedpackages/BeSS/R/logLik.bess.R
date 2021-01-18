logLik.bess=function(object,...){
  n=object$nsample
  if(object$family!="bess_gaussian"){
    deviance=-object$deviance/2
    nulldeviance=-object$nulldeviance/2
    out=c(nulldeviance, deviance)
  }else{
    deviance=-n*log(object$mse)/2
    nulldeviance=-n*log(object$nullmse)/2
    out=c(nulldeviance, deviance)
  }
  names(out)=c('nullLoglik',colnames(object$beta))
  return(out)

}

logLik.bess.one=function(object,...){
  n=object$nsample
  if(object$family!="bess_gaussian"){
    deviance=-object$deviance/2
    nulldeviance=-object$nulldeviance/2
    out=c(nulldeviance, deviance)
  }else{
    deviance=-n*log(object$mse)/2
    nulldeviance=-n*log(object$nullmse)/2
    out=c(nulldeviance, deviance)
  }
  names(out)=c('nullLoglik','Loglik')
  return(out)

}
