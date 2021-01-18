bess.cox=function(x, y, beta0, s, cox.max=20, max.steps=20, factor=NULL,
                  weights=rep(1,nrow(x)), normalize=FALSE)
{
  if(missing(beta0)) beta0=rep(0,ncol(x))
  if(is.matrix(y)!=0) y=as.matrix(y)
  if(ncol(y)!=2) stop("Please input y with two columns!")
  if(s>length(beta0))
  {stop("s is too large")}
  if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x))
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
    fit = gbess.cox(x, y, Gi, beta0=beta0, s = s,
                    max.steps = max.steps, cox.max = cox.max,
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
      mark=order(y[,1],decreasing = FALSE)
      y=y[mark,]
      x=x[mark,]
      weights=weights[mark]


      one = rep(1, n)
      #center
      meanx = drop(weights %*% x)/n
      x = scale(x, meanx, FALSE)
      #normalize
      normx = sqrt(drop(weights %*% (x^2)))
      nosignal = normx/sqrt(n) < .Machine$double.eps
      if (any(nosignal))  normx[nosignal] = (.Machine$double.eps) * sqrt(n)

      names(normx) = NULL
      x = sqrt(n)*scale(x, FALSE, normx)
    }


    ind=which(y[,2]==0)

    setA=getcox_A(x,beta,s,rep(0,p),status=ind,weights=weights)
    l=1

    A=list()
    I=list()

    A[[1]]=0
    I[[1]]=seq(p)
    A[[l+1]] = setA$A
    I[[l+1]] = setA$I

    while ((l <= max.steps))
    {

      beta[I[[l+1]]] = 0

      cox=coxph(Surv(y[,1],y[,2])~x[,A[[l+1]]],weights=weights,eps=1e-8,iter.max=cox.max)
      beta[A[[l+1]]]=cox$coefficients

      setA=getcox_A(x,beta,s,A[[l+1]],status=ind,weights=weights)

      A[[l+2]] = setA$A
      I[[l+2]] = setA$I

      if(setequal(A[[l+2]],A[[l]])|setequal(A[[l+2]],A[[l+1]])) {break}
      else{l=l+1
      gc()}
    }

    names(beta) = vn
    xbest=xs[,which(beta!=0)]
    bestmodel=coxph(Surv(y[,1],y[,2])~xbest, weights=weights, iter.max=cox.max)

    dev=-2*cox$loglik[2]
    nulldev=-2*cox$loglik[1]
    aic=dev+2*s
    bic=dev+log(n)*s
    ebic=dev+(log(n)+2*log(p))*s

    if(normalize)
    {
      beta=sqrt(n)*beta/normx
    }

    return(list(family="bess_cox",beta=beta,nsample=n,deviance=dev,bestmodel=bestmodel,
                nulldeviance=nulldev,lambda=setA$max_T^2/2,AIC=aic,BIC=bic,EBIC=ebic,max.steps=max.steps))
  }

}
