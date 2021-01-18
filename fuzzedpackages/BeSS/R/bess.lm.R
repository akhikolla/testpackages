bess.lm=function(x, y, beta0, s, max.steps=20, factor = NULL,
                 weights=rep(1,nrow(x)), normalize=FALSE)
{
  if(missing(beta0)) beta0=rep(0,ncol(x))
  if(s>length(beta0))
  {stop("s is too large")}
  if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x),"g")
  if(!is.null(factor)){
    if(!is.data.frame(x)) x = as.data.frame(x)
    x[,factor] = apply(x[,factor,drop=FALSE], 2, function(x){
      x = as.factor(x)
    })
    group = rep(1, ncol(x))
    names(group) = colnames(x)
    group[factor] = apply(x[,factor,drop=FALSE], 2,function(x) {length(unique(x))})-1
    Gi = rep(1:ncol(x), times = group)
    beta0 = rep(beta0, times = group)
    x = model.matrix(~., data = x)[,-1]
    fit = gbess.lm(x, y, Gi, beta0, s = s, max.steps = max.steps,
                   weights = weights, normalize = normalize)
    fit$factor = factor
    return(fit)
  }else{
    x = as.matrix(x)
    n = dim(x)[1]
    p = dim(x)[2]
    vn = dimnames(x)[[2]]
    one = rep(1,n)
    beta=beta0
    names(beta) = vn
    xs=x
    ys=y

    weights = weights/mean(weights)
    if(normalize)
    {
      #center
      meanx = drop(weights %*% x)/n
      x = scale(x, meanx, FALSE)
      mu = mean(y*weights)
      y = drop(y - mu)
      #normalize
      normx = sqrt(drop(weights %*% (x^2)))
      nosignal = normx/sqrt(n) < (.Machine$double.eps)
      if (any(nosignal))  normx[nosignal] = (.Machine$double.eps) * sqrt(n)

      names(normx) = NULL
      x = sqrt(n)*scale(x, FALSE, normx)
    }

    fit=bess_lm(x*sqrt(weights),y*sqrt(weights),s,max.steps,beta0)
    beta=fit$beta
    names(beta) = vn
    xbest=xs[,which(beta!=0)]
    bestmodel=lm(ys~xbest, weights = weights)

    lambda=fit$max_T^2/2
    mse=mean(weights*(y-x%*%beta)^2)
    nullmse=mean(weights*(y^2))
    aic=n*log(mse)+2*s
    bic=n*log(mse)+log(n)*s
    ebic=n*log(mse)+(log(n)+2*log(p))*s
    if(normalize)
    {
      beta=sqrt(n)*beta/normx
      coef0=mu-sum(beta*meanx)
      return(list(family="bess_gaussian",beta=beta,coef0=coef0,nsample=n,bestmodel=bestmodel,
                  lambda=lambda,mse=mse,nullmse=nullmse,AIC=aic,BIC=bic,EBIC=ebic,max.steps=max.steps))
    }else return(list(family="bess_gaussian",beta=beta,coef0=0,nsample=n,bestmodel=bestmodel,
                      lambda=lambda,mse=mse,nullmse=nullmse,AIC=aic,BIC=bic,EBIC=ebic,max.steps=max.steps))
  }
}



