deviance.bess=function(object,...){
  n=object$nsample
  if(object$family!="bess_gaussian"){
    deviance=object$deviance
    nulldeviance=object$nulldeviance
    out=c(nulldeviance, deviance)
  }else{
    deviance=n*log(object$mse)
    nulldeviance=n*log(object$nullmse)
    out=c(nulldeviance, deviance)
  }

  names(out)=c('nulldeviance',colnames(object$beta))
  return(out)

}

deviance.bess.one=function(object,...)
{
  n=object$nsample
  if(object$family!="bess_gaussian"){
    deviance=object$deviance
    nulldeviance=object$nulldeviance
    out=c(nulldeviance, deviance)
  }else{
    deviance=n*log(object$mse)
    nulldeviance=n*log(object$nullmse)
    out=c(nulldeviance, deviance)
  }

  names(out)=c('nulldeviance','deviance')
  return(out)

}
