gbess.lm = function(x, y, Gi, beta0, s, max.steps = 20,
                    weights=rep(1,nrow(x)), normalize=FALSE)
{
  if(missing(beta0)) beta0=rep(0,ncol(x))
  if(s>length(beta0))
  {stop("s is too large")}
  # initial
  n = dim(x)[1]
  p = dim(x)[2]
  vn = dimnames(x)[[2]]
  beta=beta0
  names(beta) = vn
  xs=x
  ys=y

  weights = weights/mean(weights)

  orderGi = order(Gi)
  x = x[,orderGi]
  Gi = Gi[orderGi]
  gi = unique(Gi)
  gi_index = match(gi, Gi)
  N = length(gi)
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
  x = x*sqrt(weights)
  y = y*sqrt(weights)
  PhiG = lapply(1:N, function(i){
    idx <- which(Gi==i)
    if(length(idx) == 1)
      return(-sqrt(t(x[,idx])%*%x[,idx])) else{
        return(-EigenR(t(x[,idx])%*%x[,idx]))
      }
  })
  invPhiG = lapply(PhiG, solve)
  fit = gbess_lm(X=x, y=y, G=Gi, index=gi_index, PhiG=PhiG, invPhiG=invPhiG,
                   T0=s, max_steps = max.steps, beta0 = beta0,
                   n=n, p=p, N=N)

  if(normalize)
  {
    beta=sqrt(n)*beta/normx
    coef0=mu-sum(beta*meanx)
  }
  beta = fit$beta
  beta[orderGi] = beta
  names(beta) = vn
  A = fit$A+1
  A = orderGi[A]
  B = fit$B+1
  B = orderGi[B]
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
                lambda=lambda,mse=mse,nullmse=nullmse,AIC=aic,BIC=bic,EBIC=ebic,max.steps=max.steps,
                gr_size=fit$gr_size))
  }else return(list(family="bess_gaussian",beta=beta,coef0=0,nsample=n,bestmodel=bestmodel,
                    lambda=lambda,mse=mse,nullmse=nullmse,AIC=aic,BIC=bic,EBIC=ebic,max.steps=max.steps,
                    gr_size=fit$gr_size))
}


