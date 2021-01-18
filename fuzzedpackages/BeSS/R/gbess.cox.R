gbess.cox = function(x, y, Gi, beta0, s, cox.max=20, max.steps=10,
                     weights=rep(1,nrow(x)), normalize=FALSE)
{
  if(missing(beta0)) beta0=rep(0,ncol(x))
  if(is.matrix(y)!=0) y=as.matrix(y)
  if(ncol(y)!=2) stop("Please input y with two columns!")
  if(missing(beta0)) beta0=rep(0,ncol(x))
  if(s>length(beta0))
  {stop("s is too large")}
  # initial
  n = dim(x)[1]
  p = dim(x)[2]
  vn = dimnames(x)[[2]]
  one = rep(1,n)
  names(beta0) = vn
  xs = x
  weights = weights/mean(weights)

  orderGi = order(Gi)
  x = x[,orderGi]
  Gi = Gi[orderGi]
  gi = unique(Gi)
  gi_index = match(gi, Gi)
  N = length(gi)
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

  beta = beta0
  A0 = NULL
  B = rep(0,p+1)

  for(k in 1:max.steps){
    setA = ggetcox_A(x, Gi, gi_index, s, beta, n, p, N, y[,2], weights, B)
    A = setA$A+1
    B = setA$B+1
    beta = rep(0,p)
    gr_size = setA$gr_size
    cox=coxph(Surv(y[,1],y[,2])~x[,B],weights=weights,eps=1e-8,iter.max=cox.max)
    beta[B]=cox$coefficients
    if(setequal(A,A0)==TRUE){
      break;
    }
    A0 <- A
  }
  if(normalize)
  {
    beta=sqrt(n)*beta/normx
  }
  beta[orderGi] = beta
  names(beta) = vn
  A = orderGi[A]
  B = orderGi[B]
  s=length(B)

  xbest=xs[,which(beta!=0)]
  bestmodel=coxph(Surv(y[,1],y[,2])~xbest, weights=weights, iter.max=cox.max)
  dev=-2*cox$loglik[2]
  nulldev=-2*cox$loglik[1]
  aic=dev+2*s
  bic=dev+log(n)*s
  ebic=dev+(log(n)+2*log(p))*s

  return(list(family="bess_cox",beta=beta,nsample=n,bestmodel=bestmodel,
              deviance=dev,nulldeviance=nulldev,
              lambda=setA$max_T^2/2,AIC=aic,BIC=bic,EBIC=ebic,max.steps=max.steps,
              gr_size=gr_size))
}


