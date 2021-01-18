bess.glm=function(x, y, beta0, intercept=0, s, max.steps=20,
                  glm.max=1e6, factor = NULL,
                  weights=rep(1,nrow(x)), normalize=FALSE)
{
  if(length(unique(y))!=2)  stop("Please input binary variable!")
  if(missing(beta0)) beta0=rep(0,ncol(x))
  if(s>length(beta0))
  {stop("s is too large")}
  if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x),"g")
  if(!is.null(factor)){
    if(!is.data.frame(x)) x = as.data.frame(x)
    x[,factor] = apply(x[,factor,drop=FALSE], 2, function(x){
      return(as.factor(x))
    })
    group = rep(1, ncol(x))
    names(group) = colnames(x)
    group[factor] = apply(x[,factor,drop=FALSE], 2,function(x) {length(unique(x))})-1
    Gi = rep(1:ncol(x), times = group)
    beta0 = rep(beta0, times = group)
    x = model.matrix(~., data = x)[,-1]
    fit = gbess.glm(x, y, Gi, beta0=beta0, intercept=intercept, s = s,
                    max.steps = max.steps, glm.max=glm.max,
                    weights = weights, normalize = normalize)
    fit$factor = factor
    return(fit)
  }else{
    x = as.matrix(x)
    n = dim(x)[1]
    p = dim(x)[2]
    vn = dimnames(x)[[2]]
    one = rep(1,n)
    beta = beta0
    names(beta) = vn
    xs = x

    weights = weights/mean(weights)

    if(normalize)
    {
      meanx = drop(weights %*% x)/n
      x = scale(x, meanx, FALSE)

      normx = sqrt(drop(weights %*% (x^2)))
      nosignal = normx/sqrt(n) < .Machine$double.eps
      if (any(nosignal))  normx[nosignal] = (.Machine$double.eps) * sqrt(n)

      names(normx) = NULL
      x = sqrt(n)*scale(x, FALSE, normx)
    }
    setA=get_A(x,y,beta,intercept,s,rep(0,p),weights)
    pr=setA$p
    l=1

    A=list()
    I=list()

    A[[1]]=0
    I[[1]]=seq(p)
    A[[l+1]] = setA$A
    I[[l+1]] = setA$I
    S=1:nrow(x)

    while (l <= max.steps)
    {
      beta[I[[l+1]]] = 0

      if(s>=2)
      {
        logit=glmnet(x[,A[[l+1]]],y,family="binomial",lambda = 0,maxit=glm.max,weights=weights)
        beta[A[[l+1]]]=logit$beta
        coef0=logit$a0
      }else{
        logit=glm(y~x[,A[[l+1]]],family="binomial",weights=weights)
        beta[A[[l+1]]]=logit$coefficients[-1]
        coef0=logit$coefficients[1]
      }

      setA=get_A(x,y,beta,coef0,s,A[[l+1]],weights)
      pr=setA$p

      A[[l+2]] = setA$A
      I[[l+2]] = setA$I

      if(setequal(A[[l+2]],A[[l]])|setequal(A[[l+2]],A[[l+1]])) {break}
      else{l=l+1
      gc()}
    }
    #dev=logit$deviance
    names(beta) = vn
    xbest=xs[,which(beta!=0)]
    bestmodel=glm(y~xbest, family=binomial, weights=weights)

    dev=-2*sum((weights*((y*log(pr) + (1-y)*log(1-pr))))[which(pr>1e-20&pr<1-1e-20)])
    nulldev=-2*sum(weights*(y*log(0.5) + (1-y)*log(0.5)))
    aic=dev+2*s
    bic=dev+log(n)*s
    ebic=dev+(log(n)+2*log(p))*s

    if(normalize)
    {
      beta=sqrt(n)*beta/normx
      coef0=coef0-sum(beta*meanx)
    }

    return(list(family="bess_binomial",beta=beta,coef0=coef0,nsample=n,bestmodel=bestmodel,
                deviance=dev,nulldeviance=nulldev,
                lambda=setA$max_T^2/2,p=pr,AIC=aic,BIC=bic,EBIC=ebic,max.steps=max.steps))
  }
}





