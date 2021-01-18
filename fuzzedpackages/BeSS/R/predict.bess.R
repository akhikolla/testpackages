predict.bess=function(object, newdata, type = c("ALL", "opt", "AIC", "BIC", "EBIC"),...)
{
  type <- match.arg(type)
  method = object$method
  if(method == "gsection"&(type%in%c("AIC","BIC","EBIC"))) stop("method gsection shouldn't match type AIC, BIC or EBIC!")
  if(method == "sequential"&type=="opt") stop("method sequential shouldn't match type opt!")
  if(!is.null(object$factor)){
    factor = c(object$factor)
    if(!is.data.frame(newdata)) newdata = as.data.frame(newdata)
    newdata[,factor] = apply(newdata[,factor,drop=FALSE], 2, function(x){
      return(as.factor(x))
    })
    newdata = model.matrix(~., data = newdata)[,-1]
  }
  if(is.null(colnames(newdata))) {
    newx = as.matrix(newdata)
  }else{
    vn = rownames(object$beta)
    if(any(is.na(match(vn, colnames(newdata))))) stop("names of newdata don't match training data!")
    newx = as.matrix(newdata[,vn])
  }
  if(object$family == "bess_gaussian")
  {
    betas = object$beta
    coef0 = object$coef0
    y = t(newx%*%betas)+coef0
    if(type == "ALL"){
      return(y)
    }
    if(type == "opt"){
      return(y[nrow(y),,drop = TRUE])
    }
    return(y[which.min(object[[type]]),,drop = TRUE])
  }
  if(object$family == "bess_binomial")
  {
    betas = object$beta
    coef = object$coef0
    class = matrix(0,ncol(betas),nrow(newx))
    for(i in 1:ncol(betas))
    {
      class[i,] = as.numeric(exp(newx%*%betas[,i]+coef[i])/(1+exp(newx%*%betas[,i]+coef[i]))>0.5)
      class[i,][is.na(class[i,])] = 1
      if(!is.null(object$y_names))
      {
       class[which(class == 0,arr.ind = T)] = object$y_names[1]
       class[which(class == 1,arr.ind = T)] = object$y_names[2]
      }
    }
    if(type == "ALL"){
      return(class)
    }
    if(type == "opt"){
      return(class[nrow(class),,drop = TRUE])
    }
    return(class[which.min(object[[type]]),,drop = TRUE])
  }
  if(object$family=="bess_cox")
  {
    betas = object$beta

    betax = newx%*%betas
    if(type == "ALL"){
      return(t(betax))
    }
    if(type == "opt"){
      return(t(betax)[ncol(betax),,drop = TRUE])
    }
    return(betax[,which.min(object[[type]]),drop = TRUE])
  }

}



predict.bess.one=function(object,newdata, ...)
{
  if(!is.null(object$factor)){
    factor = c(object$factor)
    if(!is.data.frame(newdata)) newdata = as.data.frame(newdata)
    newdata[,factor] = apply(newdata[,factor,drop=FALSE], 2, function(x){
      return(as.factor(x))
    })
    newdata = model.matrix(~., data = newdata)[,-1]
  }
  if(is.null(colnames(newdata))) {
    newx = as.matrix(newdata)
  }else{
    vn = names(object$beta)
    if(any(is.na(match(vn, colnames(newdata))))) stop("names of newdata don't match training data!")
    newx = as.matrix(newdata[,vn])
  }
  if(object$family=="bess_gaussian")
  {
    betas = object$beta
    coef0 = object$coef0

    y = drop(newx %*% betas)+coef0
    return(y)
  }
  if(object$family == "bess_binomial")
  {
    betas = object$beta
    coef = object$coef0

    class = as.numeric(exp(newx%*%betas+coef)/(1+exp(newx%*%betas+coef))>0.5)
    class[is.na(class)] = 1
    if(!is.null(object$y_names))
    {
      class[which(class == 0,arr.ind = T)] = object$y_names[1]
      class[which(class == 1,arr.ind = T)] = object$y_names[2]
    }

    return(class)
  }
  if(object$family == "bess_cox")
  {
    betas = object$beta

    betax = newx%*%betas;
    return(drop(betax))
  }

}

