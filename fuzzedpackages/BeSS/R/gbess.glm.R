gbess.glm = function(x, y, Gi, beta0, intercept=0, s, max.steps = 10, glm.max=1e6,
                     weights=rep(1,nrow(x)), normalize=FALSE)
{
  if(length(unique(y))!=2)  stop("Please input binary variable!")

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
    meanx = drop(weights %*% x)/n
    x = scale(x, meanx, FALSE)

    normx = sqrt(drop(weights %*% (x^2)))
    nosignal = normx/sqrt(n) < .Machine$double.eps
    if (any(nosignal))  normx[nosignal] = (.Machine$double.eps) * sqrt(n)

    names(normx) = NULL
    x = sqrt(n)*scale(x, FALSE, normx)
  }

  beta = beta0
  coef0 = intercept
  A0 = NULL
  B = rep(0,p+1)

  for(k in 1:max.steps){
    setA = gget_A(x, y, Gi, gi_index, s, beta, coef0, n, p, N, weights, B)
    A = setA$A+1
    B = setA$B+1
    beta = rep(0,p)
    gr_size = setA$gr_size
    if(length(B)>=2)
    {
      logit=glmnet(x[,B],y,family="binomial",lambda = 0,maxit=glm.max, weights = weights)
      beta[B]=logit$beta
      coef0=logit$a0
    }else{
      logit=glm(y~x[,B],family="binomial", weights = weights)
      beta[B]=logit$coefficients[-1]
      coef0=logit$coefficients[1]
    }
    if(setequal(A,A0)==TRUE){
      break;
    }
    A0 <- A
  }
  if(normalize)
  {
    beta=sqrt(n)*beta/normx
    coef0=coef0-sum(beta*meanx)
  }
  beta[orderGi] = beta
  names(beta) = vn
  A = orderGi[A]
  B = orderGi[B]
  s=length(B)
  eta = x%*%beta
  pr = exp(eta)/(1+exp(eta))

  xbest=xs[,which(beta!=0)]
  bestmodel=glm(y~xbest, family="binomial", weights = weights)
  dev=-2*sum((weights*((y*log(pr) + (1-y)*log(1-pr))))[which(pr>1e-20&pr<1-1e-20)])
  nulldev=-2*sum(weights*(y*log(0.5) + (1-y)*log(0.5)))
  aic=dev+2*s
  bic=dev+log(n)*s
  ebic=dev+(log(n)+2*log(p))*s

  return(list(family="bess_binomial",beta=beta,coef0=coef0,nsample=n,bestmodel=bestmodel,
              deviance=dev,nulldeviance=nulldev,
              lambda=setA$max_T^2/2,p=p,AIC=aic,BIC=bic,EBIC=ebic,max.steps=max.steps,
              gr_size=gr_size))
}


