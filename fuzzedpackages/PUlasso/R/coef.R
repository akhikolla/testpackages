#' @export
#' @method coef PUfit
coef.PUfit<-function(object,lambda=NULL,std.scale=F,...){
  fit = object
  if(std.scale){
    fitcoef = fit$std_coef
  }else{
    fitcoef = fit$coef
  }
  
  
  if(is.null(lambda)){
    return(fitcoef)
  }else{
    nlambda = length(lambda)
    coef = matrix(0, nrow = nrow(fitcoef),ncol=nlambda)
    rownames(coef) = rownames(fitcoef)
    for(i in 1:nlambda){
      lidx<-which(fit$lambda==lambda[i])
      if(length(lidx)==0){
        if(max(fit$lambda)<lambda[i]){
          coef[,i] = fitcoef[,1]
        }else if(min(fit$lambda)>lambda[i]){
          coef[,i] = fitcoef[,length(fit$lambda)]
        }else{
          lidx1<-max(which(fit$lambda>lambda[i]))
          lidx2<-min(which(fit$lambda<=lambda[i]))
          t = (lambda[i] - fit$lambda[lidx2])/(fit$lambda[lidx1] - fit$lambda[lidx2])
          coef[,i] = fitcoef[,lidx2]+t*(fitcoef[,lidx1]-fitcoef[,lidx2])
        }
      }else{
        coef[,i] = fitcoef[,lidx]
      }
    }
  }
  return(coef)
}

#' @export
#' @method coef cvPUfit
coef.cvPUfit<-function(object,lambda=NULL,std.scale=F,...){
  fit = object$PUfit
  
  if(std.scale){
    fitcoef = fit$std_coef
  }else{
    fitcoef = fit$coef
  }
  
  
  if(is.null(lambda)){
    return(fitcoef)
  }else{
    nlambda = length(lambda)
    coef = matrix(0, nrow = nrow(fitcoef),ncol=nlambda)
    rownames(coef) = rownames(fitcoef)
    for(i in 1:nlambda){
      lidx<-which(fit$lambda==lambda[i])
      if(length(lidx)==0){
        if(max(fit$lambda)<lambda[i]){
          coef[,i] = fitcoef[,1]
        }else if(min(fit$lambda)>lambda[i]){
          coef[,i] = fitcoef[,length(fit$lambda)]
        }else{
          lidx1<-max(which(fit$lambda>lambda[i]))
          lidx2<-min(which(fit$lambda<=lambda[i]))
          t = (lambda[i] - fit$lambda[lidx2])/(fit$lambda[lidx1] - fit$lambda[lidx2])
          coef[,i] = fitcoef[,lidx2]+t*(fitcoef[,lidx1]-fitcoef[,lidx2])
        }
      }else{
        coef[,i] = fitcoef[,lidx]
      }
    }
  }
  return(coef)
}
  