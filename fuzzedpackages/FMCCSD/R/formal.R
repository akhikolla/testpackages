.quadrature.rules <- function( recurrences, inner.products )
{
  
  np1 <- nrow( recurrences )
  n <- np1 - 1
  rules <- as.list( rep( NULL, n ) )
  monic.recurrences <- orthopolynom::monic.polynomial.recurrences( recurrences )
  matrices <- orthopolynom::jacobi.matrices( monic.recurrences )
  matrix.eigens <- lapply( matrices, eigen )
  roots <- orthopolynom::polynomial.roots( monic.recurrences )
  h.0 <- inner.products[1]
  for ( k in 1:n ) {
    values <- matrix.eigens[[k]]$values
    vectors <- matrix.eigens[[k]]$vectors
    x <- values
    w <- rep( 0, k )
    for ( j in 1:k ) {
      v.j <- vectors[1,j]
      w[j] <- h.0 * v.j * v.j
    }
    rule <- data.frame( cbind( x, w ) )
    names( rule ) <- c( "x", "w" )
    rules[[k]] <- rule
  }
  return( rules )
}


.ghrules=function(n,normalized=FALSE){
  r=orthopolynom::hermite.h.recurrences(n,normalized)
  ip=orthopolynom::hermite.h.inner.products(n)
  return(.quadrature.rules(r,ip))
}

.Data.trans=function(Rawdata,n_subject.cov,n_tooth.cov){
  n=length(unique(Rawdata[,1]))
  ni=sort(unique(Rawdata[,1]))
  mi=rep(0,n)
  for(i in 1:n){
    mi[i]=sum(Rawdata[,1]==ni[i])
  }
  CSTime=list()
  length(CSTime)=n
  X=list()
  length(X)=n
  Delta=list()
  length(Delta)=n
  for (i in 1:n) {
    CSTime[[i]]=Rawdata[,2][Rawdata[,1]==ni[i]]
    X[[i]]=as.matrix(Rawdata[Rawdata[,1]==ni[i],(3+n_subject.cov):(dim(Rawdata)[2]-1)])
    Delta[[i]]=Rawdata[,dim(Rawdata)[2]][Rawdata[,1]==ni[i]]
  }
  Z=matrix(0, nrow = n, ncol = n_subject.cov)
  
  for (i in 1:n) {
    
    Z[i,]=as.numeric(Rawdata[sum(mi[0:i]),3:(2+n_subject.cov)])
  }
  Z=as.matrix(Z)
  
  finalresult=list()
  length(finalresult)=6
  finalresult[[1]]=n
  finalresult[[2]]=mi
  finalresult[[3]]=CSTime
  finalresult[[4]]=Delta
  finalresult[[5]]=X
  finalresult[[6]]=Z
  return(finalresult)
}
.Datascale=function(Rawdata,n_subject.cov,n_tooth.cov){
  covdata=Rawdata[,3:(2+n_subject.cov+n_tooth.cov)]
  covdim=dim(covdata)[2]
  
  for(i in 1:covdim){
    if(is.numeric(covdata[,i])){
      covdata[,i]=scale(covdata[,i])[,1]
    }
    else{
      covdata[,i]=as.numeric(covdata[,i])-1
    }
  }
  finalresult=cbind(Rawdata[,1:2],covdata,Rawdata[,dim(Rawdata)[2]])
  return(finalresult)
}

CSDfit=function(Rawdata,n_subject.raw,n_within.raw,r,n_quad=30,lambda=0,tolerance=0.5,
                   knots.num=2,degree=2,scale.numr=TRUE,clustering=TRUE){
  myrules=.ghrules(n_quad,normalized=FALSE)
  myrules=as.matrix(myrules[[n_quad]])
  if(scale.numr==TRUE){
    Rawdata=.Datascale(Rawdata,n_subject.raw,n_within.raw)
    n_tooth.cov=n_within.raw
    n_subject.cov=n_subject.raw
    totaldata=.Data.trans(Rawdata,n_subject.cov,n_tooth.cov)
  }
  else{
    totaldata=.Data.trans(Rawdata,n_subject.raw,n_within.raw)
    n_tooth.cov=n_within.raw
    n_subject.cov=n_subject.raw
  }
  n=totaldata[[1]]
  ni=totaldata[[2]]
  C=totaldata[[3]]
  Delta=totaldata[[4]]
  X=totaldata[[5]]
  Z=totaldata[[6]]
  minCSTime=min(Rawdata[,2])
  maxCSTime=max(Rawdata[,2])
  blC <- list()
  length(blC) <- n
  knots <- seq(0,1  , length.out = (knots.num + 2))
  knots=knots[3:length(knots)-1]
  for (i in 1:n) {
    blC[[i]]=t(splines2::ibs((C[[i]]-(minCSTime-0.1))/(maxCSTime-minCSTime+0.2),knots = knots,degree=degree
                             ,Boundary.knots = c(0,1),intercept = TRUE))
  }
  if(clustering){
    coefnum=dim(blC[[1]])[1]
    D=matrix(0, nrow = coefnum-2, ncol = coefnum)
    for(i in 1:(coefnum-2)){
      D[i,i:(i+2)]=c(1,-2,1)
    }
    R=(t(D)%*%D)
    betadim=dim(X[[1]])[2]
    gammadim=dim(Z)[2]
    lastpar=rep(0,betadim+gammadim+1+knots.num+degree+1)
    difference=1
    while (difference>tolerance) {
      output=optim(par=rep(0,length(lastpar)),fn=targetfunc,lastpar=lastpar,rules=myrules,C=C,Delta=Delta,X=X,Z=Z,n=n,ni=ni,r=r,blC=blC,betadim=betadim,
                    gammadim=gammadim,R=R,lambda=lambda,method = "BFGS")
      outputpar=output$par
      difference=sum(abs((-outputpar[1:(betadim+gammadim+1)]+lastpar[1:(betadim+gammadim+1)])/(lastpar[1:(betadim+gammadim+1)])))
      lastpar=outputpar
    }
    
    hessian=numDeriv::hessian(func=testquadrature1current,x=outputpar,rules=myrules,C=C,Delta=Delta,
                              X=X,Z=as.matrix(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim,gammadim=gammadim,R=R,lambda=lambda)
    var=diag(solve(hessian,tol=1e-40))
    
    parest=outputpar
    pardim=length(parest)
    
    psi=parest[(betadim+gammadim+2):length(parest)]
    Rbig=matrix(0, nrow = length(parest), ncol = length(parest))
    Rsecderiv=numDeriv::hessian(penaltyterm,x=psi,lambda=lambda,R=R)
    Rbig[(betadim+gammadim+2):length(parest),(betadim+gammadim+2):length(parest)]=Rsecderiv
    df=sum(diag((hessian-Rbig)%*%solve(hessian,tol=1e-50)))
    minusloglikelihood=testquadrature1current(parest,rules=myrules,C=C,Delta=Delta,
                                              X=X,Z=as.matrix(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim
                                              ,gammadim=gammadim,R=R,lambda=lambda)
    AIC=df+minusloglikelihood
    
    reg.est=parest[1:(n_subject.cov+n_tooth.cov+1)]
    reg.est[length(reg.est)]=exp(reg.est[length(reg.est)])
    reg.se=sqrt(var[1:(n_subject.cov+n_tooth.cov+1)])
    reg.se[length(reg.se)]=reg.est[length(reg.est)]*reg.se[length(reg.se)]
    z.stat=abs(reg.est/reg.se)
    p.est=2*(1-pnorm(z.stat))
    resultmat=matrix(0, nrow = n_subject.cov+n_tooth.cov+1, ncol = 6)
    resultmat[,1]=reg.est
    resultmat[,2]=reg.se
    resultmat[,3]=z.stat
    resultmat[,4]=p.est
    resultmat[,5]=resultmat[,1]-1.96*resultmat[,2]
    resultmat[,6]=resultmat[,1]+1.96*resultmat[,2]
    resultmat[length(reg.est),5]=exp(parest[length(reg.est)]-1.96*sqrt(var[length(reg.est)]))
    resultmat[length(reg.est),6]=exp(parest[length(reg.est)]+1.96*sqrt(var[length(reg.est)]))
    resultmat=round(resultmat,2)
    resultmat=data.frame(resultmat)
    colnames(resultmat)=c("par.est","SE","Z","p-value","CI_lower","CI_upper")
    rownames(resultmat)[dim(resultmat)[1]]="Random_effect"

    for(k in 1:n_tooth.cov){
      rownames(resultmat)[k]=paste("Surv",colnames(Rawdata)[2+n_subject.cov+k],sep = "_")
    }
    for (k in 1:n_subject.cov) {
      rownames(resultmat)[n_tooth.cov+k]=paste("Surv",colnames(Rawdata)[2+k],sep = "_")
    }
    coefs=round(exp(parest[(n_subject.cov+n_tooth.cov+2):length(parest)]),2)
  }
  else{
    coefnum=dim(blC[[1]])[1]
    D=matrix(0, nrow = coefnum-2, ncol = coefnum)
    for(i in 1:(coefnum-2)){
      D[i,i:(i+2)]=c(1,-2,1)
    }
    R=(t(D)%*%D)
    betadim=dim(X[[1]])[2]
    gammadim=dim(Z)[2]
    lastpar=rep(0,betadim+gammadim+knots.num+degree+1)
    difference=1
    while (difference>tolerance) {
      output=optim(par=rep(0,length(lastpar)),fn=targetfuncfrailty,lastpar=lastpar,rules=myrules,C=C,Delta=Delta,X=X,Z=Z,n=n,ni=ni,r=r,blC=blC,betadim=betadim,
                            gammadim=gammadim,R=R,lambda=lambda,method = "BFGS")
      outputpar=output$par
      difference=sum(abs((-outputpar[1:(betadim+gammadim)]+lastpar[1:(betadim+gammadim)])/(lastpar[1:(betadim+gammadim)])))
      lastpar=outputpar
    }
    hessian=numDeriv::hessian(func=testquadrature1currentnofrailty,x=outputpar,rules=myrules,C=C,Delta=Delta,
                              X=X,Z=as.matrix(Z),n=n,ni=ni,r=r,blC=blC,betadim=1,gammadim=1,Cauchyindex=FALSE,Cauchyscale=2.5)
    var=diag(solve(hessian,tol=1e-40))
    
    parest=outputpar
    pardim=length(parest)
    
    psi=parest[(betadim+gammadim+1):length(parest)]
    Rbig=matrix(0, nrow = length(parest), ncol = length(parest))
    Rsecderiv=numDeriv::hessian(penaltyterm,x=psi,lambda=lambda,R=R)
    Rbig[(betadim+gammadim+1):length(parest),(betadim+gammadim+1):length(parest)]=Rsecderiv
    df=sum(diag((hessian-Rbig)%*%solve(hessian,tol=1e-50)))
    minusloglikelihood=testquadrature1currentnofrailty(parest,rules=myrules,C=C,Delta=Delta,
                                              X=X,Z=as.matrix(Z),n=n,ni=ni,r=r,blC=blC,betadim=betadim
                                              ,gammadim=gammadim,R=R,lambda=lambda)
    
    reg.est=parest[1:(n_subject.cov+n_tooth.cov)]
    reg.se=sqrt(var[1:(n_subject.cov+n_tooth.cov)])
    z.stat=abs(reg.est/reg.se)
    p.est=2*(1-pnorm(z.stat))
    resultmat=matrix(0, nrow = n_subject.cov+n_tooth.cov+1, ncol = 6)
    resultmat[,1]=reg.est
    resultmat[,2]=reg.se
    resultmat[,3]=z.stat
    resultmat[,4]=p.est
    resultmat[,5]=resultmat[,1]-1.96*resultmat[,2]
    resultmat[,6]=resultmat[,1]+1.96*resultmat[,2]
    resultmat=round(resultmat,2)
    resultmat=data.frame(resultmat)
    colnames(resultmat)=c("par.est","SE","Z","p-value","CI_lower","CI_upper")

    for(k in 1:n_tooth.cov){
      rownames(resultmat)[k]=paste("Surv",colnames(Rawdata)[2+n_subject.cov+k],sep = "_")
    }
    for (k in 1:n_subject.cov) {
      rownames(resultmat)[n_tooth.cov+k]=paste("Surv",colnames(Rawdata)[2+k],sep = "_")
    }
    coefs=round(exp(parest[(n_subject.cov+n_tooth.cov+1):length(parest)]),2)
  }
  
  
  return(list(parameter.est=resultmat,log_likelihood=-minusloglikelihood,AICvalue=AIC,
              coefs=coefs))
}