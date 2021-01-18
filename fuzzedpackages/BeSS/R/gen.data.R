gen.data=function(n, p, family, K, rho = 0, sigma = 1, beta = NULL,
                  censoring = TRUE, c = 1, scal)
{
  one=rep(1,n)
  zero=rep(0,n)
  X=rnorm(n*p)
  X=matrix(X,n,p)
  X = scale(X, TRUE, FALSE)
  normX = sqrt(drop(one %*% (X^2)))
  X = sqrt(n)*scale(X, FALSE, normX)
  gc()


  x=X+rho*(cbind(zero,X[,1:(p-2)],zero)+cbind(zero,X[,3:p],zero))
  colnames(x)=paste0('X',1:ncol(x))


  rm(X)
  gc()

  nonzero=sample(1:p,K)
  Tbeta=rep(0,p)

  if(family=="gaussian")
  {
    m=5*sqrt(2*log(p)/n)
    M=100*m
    if(is.null(beta)) Tbeta[nonzero]=runif(K,m,M) else Tbeta=beta
    y=drop(x %*% Tbeta+rnorm(n,0,sigma^2))

    return(list(x=x,y=y,Tbeta=Tbeta))
  }

  if(family=="binomial")
  {
    m=5*sigma*sqrt(2*log(p)/n)
    if(is.null(beta)) Tbeta[nonzero]=runif(K,2*m,10*m) else Tbeta=beta

    ex=exp(drop(x %*% Tbeta))
    logit=ex/(1+ex)
    y=rbinom(n=n,size=1,prob=logit)

    return(list(x=x,y=y,Tbeta=Tbeta))
  }

  if(family=="cox")
  {
    m=5*sigma*sqrt(2*log(p)/n)
    if(is.null(beta)) Tbeta[nonzero]=runif(K,2*m,10*m) else Tbeta=beta


    time = (-log(runif(n))/drop(exp(x%*%Tbeta)))^(1/scal)
    if (censoring) {
      ctime = c*runif(n)
      status = (time < ctime) * 1
      censoringrate = 1 - sum(status)/n
      cat("censoring rate:", censoringrate, "\n")
      time = pmin(time, ctime)
    }else {
      status = rep(1, times = n)
      cat("no censoring", "\n")
    }

    return(list(x=x,y=cbind(time,status),Tbeta=Tbeta))
  }

}
