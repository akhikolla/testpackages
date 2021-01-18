FASTmrEMMA<-function(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,Likelihood,CLO){

if(Likelihood=="REML"){
  flagREMLE<-1
}else if(Likelihood=="ML"){
  flagREMLE<-0
}

inputform<-Genformat

if(is.null(kk)){
  
  emma.kinship <- function(snps, method="additive", use="all") {
    n0 <- sum(snps==0,na.rm=TRUE)
    nh <- sum(snps==0.5,na.rm=TRUE)
    n1 <- sum(snps==1,na.rm=TRUE)
    nNA <- sum(is.na(snps))
    #stopifnot(n0+nh+n1+nNA == length(snps))
    if ( method == "dominant" ) {
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
      snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
    }
    else if ( method == "recessive" ) {
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
      snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
    }
    else if ( ( method == "additive" ) && ( nh > 0 ) ) {
      dsnps <- snps
      rsnps <- snps
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
      #los<-intersect(which(!is.na(snps)),which(snps==0.5))
      dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
      rm(flags)
      gc()
      flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
      rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
      rm(flags,snps)
      gc()
      snps <- rbind(dsnps,rsnps)
      rm(dsnps,rsnps)
      gc()
    }
    if ( use == "all" ) {
      mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
      #losna<-which(is.na(snps))
      snps[is.na(snps)] <- mafs[is.na(snps)]
      rm(mafs)
      gc()
    }
    else if ( use == "complete.obs" ) {
      snps <- snps[rowSums(is.na(snps))==0,]
    }
    n <- ncol(snps)
    #K<-(t(snps)%*%snps+t(1-snps)%*%(1-snps))/nrow(snps)
    K<-(mrMLM::multiplication_speed(t(snps),snps)+mrMLM::multiplication_speed(t(1-snps),(1-snps)))/nrow(snps)
    diag(K) <- 1
    return(K)
  }
  
  
  if(is.null(gen)==TRUE)
  {
    warning("Please input correct genotype dataset !")
  }else{
    snp8<-gen[,3:ncol(gen)]
    kk<-emma.kinship(snp8) 
    rm(snp8)
    gc()
  }
}

if(is.null(psmatrix)){
  flagps<-1
}else{
  
  flagps<-0
}


if(is.null(svpal)==TRUE||is.null(svmlod)==TRUE){
  warning("Please set parameter!")
}

if((svpal<0)||(svpal>1))
{
  warning("Please input critical P-value between 0 and 1!")
}

if(svmlod<0)
{
  warning("Please input critical LOD score: >0!")
}
if(exists("gen")==FALSE)
{
  warning("Please input correct genotype dataset !")
}
if(exists("phe")==FALSE)
{
  warning("Please input correct phenotype dataset !")
}
if(exists("kk")==FALSE)
{
  warning("Please input correct kinship (K) dataset !")
}
if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(ncol(gen)!=(nrow(phe)+2)))
{
  warning("Sample sizes between genotypic and phenotypic datasets do not equal !")
}

if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2)))&&(svpal>=0)&&(svpal<=1)&&(svmlod>=0))
{

parmsShow=NULL
wan=NULL
parms=NULL
ress1=NULL
mannewp=NULL

multinormal<-function(y,mean,sigma)
{
  pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
  return (pdf_value)
}

ebayes_EM<-function(x,z,y)
{
  n<-nrow(z);k<-ncol(z)
  
  if(abs(min(eigen(crossprod(x,x))$values))<1e-6){
    b<-solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y)
  }else{
    b<-solve(crossprod(x,x))%*%(crossprod(x,y))
  }
  
  v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
  u<-matrix(rep(0,k),k,1)
  v<-matrix(rep(0,k),k,1)
  s<-matrix(rep(0,k),k,1)
  for(i in 1:k)
  {
    zz<-z[,i]
    s[i]<-((crossprod(zz,zz)+1e-100)^(-1))*v0
    u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
    v[i]<-u[i]^2+s[i]
  }
  
  vv<-matrix(rep(0,n*n),n,n);
  for(i in 1:k)
  {
    zz<-z[,i]
    vv=vv+tcrossprod(zz,zz)*v[i]
  }
  vv<-vv+diag(n)*v0
  
  iter<-0;err<-1000;iter_max<-500;err_max<-1e-8
  tau<-0;omega<-0
  while((iter<iter_max)&&(err>err_max))
  {
    iter<-iter+1
    v01<-v0
    v1<-v
    b1<-b
    vi<-solve(vv)
    xtv<-crossprod(x,vi)
    
    if(ncol(x)==1)
    {
      b<-((xtv%*%x)^(-1))*(xtv%*%y)
    }else{
      if(abs(min(eigen(xtv%*%x)$values))<1e-6){
        b<-solve((xtv%*%x)+diag(ncol(x))*1e-8)%*%(xtv%*%y)
      }else{
        b<-solve(xtv%*%x)%*%(xtv%*%y)
      }
    }
    r<-y-x%*%b
    ss<-matrix(rep(0,n),n,1)
    for(i in 1:k)
    {
      zz<-z[,i]
      zztvi<-crossprod(zz,vi)
      u[i]<-v[i]*zztvi%*%r
      s[i]<-v[i]*(1-zztvi%*%zz*v[i])
      v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
      ss<-ss+zz*u[i]
    }
    v0<-as.numeric(crossprod(r,(r-ss))/n)
    
    vv<-matrix(rep(0,n*n),n,n)
    for(i in 1:k)
    {
      zz<-z[,i]
      vv<-vv+tcrossprod(zz,zz)*v[i]
    }
    vv<-vv+diag(n)*v0
    
    err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
    beta<-t(b)
    sigma2<-v0
  }
  
  wang<-matrix(rep(0,k),k,1)
  for (i in 1:k){
    stderr<-sqrt(s[i]+1e-20)
    t<-abs(u[i])/stderr
    f<-t*t
    p<-pchisq(f,1,lower.tail = F)
    wang[i]<-p
  }
  
  return(list(u=u,sigma2=sigma2,wang=wang))
}

likelihood<-function(xxn,xxx,yn,bbo)
{
  nq<-ncol(xxx)
  ns<-nrow(yn)
  at1<-0
  
  if(is.null(bbo)==TRUE){
    ww1<-1:ncol(xxx)
    ww1<-as.matrix(ww1)
  }else{
    ww1<-as.matrix(which(abs(bbo)>1e-5))
  }
  at1<-dim(ww1)[1]
  lod<-matrix(rep(0,nq),nq,1)
  if(at1>0.5)
    ad<-cbind(xxn,xxx[,ww1])
  else
    ad<-xxn
  if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6)
    bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
  else
    bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
  vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns);
  ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
  sub<-1:ncol(ad);
  if(at1>0.5)
  {
    for(i in 1:at1)
    {
      ij<-which(sub!=sub[i+ncol(xxn)])
      ad1<-ad[,ij]
      if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6)
        bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
      else
        bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn) 
      vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
      ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
      lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
    }
  }
  return (lod)
}

emma.eigen.L <- function(Z,K,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.L.wo.Z(K))
  }
  else {
    return(emma.eigen.L.w.Z(Z,K,complete))
  }
}
#likelihood
emma.eigen.L.wo.Z <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}
#likelihood
emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
  return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
}

#restricted likelihood
emma.eigen.R <- function(Z,K,X,complete=TRUE) {
  if ( ncol(X) == 0 ) {
    return(emma.eigen.L(Z,K))
  }
  else if ( is.null(Z) ) {
    return(emma.eigen.R.wo.Z(K,X))
  }
  else {
    return(emma.eigen.R.w.Z(Z,K,X,complete))
  }
}
#restricted likelihood
emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}
#restricted likelihood
emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)
  
  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)    
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1:(t-q)],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                           complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
}

emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(delta*lambda+1))))-sum(log(delta*xi+1))) )  
}

emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*xi.1+1)) ))
  
}

emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- delta*lambda+1
  return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(xi/(delta*xi+1))) )
}

emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- delta*lambda+1
  return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(xi.1/(delta*xi.1+1))) )
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(delta*lambda+1))))-sum(log(delta*lambda+1))) )
}

emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- delta*lambda+1
  return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(lambda/ldelta)) )
}

emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
  t <- t1
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- delta*lambda+1
  return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
}

emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                     esp=1e-10, eig.L = NULL, eig.R = NULL)
{
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(ML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z(K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,n-q,m)    
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
    Xis.1<-matrix(eig.L$values,n,m)
    Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1
    Etasq <- matrix(etas*etas,n-q,m)
    dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
      }
    }
    
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z(Z,K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1<-matrix(eig.R$values,t-q,m)
    Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
    
    Xis.1<-matrix(eig.L$values,t,m)
    Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    
    dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
      }
    }
    
  }
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  optLL=replaceNaN(optLL)  #20160728
  
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
  }
  maxvg <- maxve*maxdelta
  
  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
  
}
emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,n-q,m)
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
    Etasq <- matrix(etas*etas,n-q,m)
    
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
      }
    }
    
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas.1 <- matrix(eig.R$values,t-q,m) 
    Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
      }
    }
    
  }  
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  optLL=replaceNaN(optLL)  
  maxLL <- max(optLL)
  
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
  }
  maxvg <- maxve*maxdelta
  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
}
################################################
#likelihood
FASTmrEMMA.delta.ML.LL.c<-function(logdelta,X,M,M.y,yMy,n){
  #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
  #n<-dim(M)[1]
  delta <-  exp(logdelta)
  ci<-as.numeric(crossprod(X))
  delta1<-as.numeric(t(X)%*%M%*%X)
  xMy<-as.numeric(crossprod(X,M.y))
  return(0.5*(n*((log(n/(2*pi))-log(as.numeric(yMy)-delta*(xMy)^2/(1+delta*delta1)))-1)-log(delta*ci+1)))
}
#dML

FASTmrEMMA.delta.ML.dLL.c<-function(logdelta,X,M,M.y,yMy,n){
  #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
  #n<-dim(M)[1]
  delta <-  exp(logdelta)
  ci<-as.numeric(crossprod(X))
  delta1<-as.numeric(t(X)%*%M%*%X)
  xMy<-as.numeric(crossprod(X,M.y))
  return(-0.5*ci/(1+delta*ci)+0.5*n/((1+delta*delta1)*(as.numeric(yMy)*(1+delta*delta1)/(xMy^2)-delta)))
}
#restrict likelihood 20190902
FASTmrEMMA.delta.REML.LL.c<-function(logdelta,X,M,M.y,yMy,v){
  #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
  #v<-n-1
  delta <-  exp(logdelta)
  #ci<-crossprod(X)
  delta1<-as.numeric(t(X)%*%M%*%X)
  xMy<-as.numeric(crossprod(X,M.y))
  return(0.5*(v*((log(v/(2*pi))-log(as.numeric(yMy)-delta*(xMy)^2/(1+delta*delta1)))-1)-log(delta*delta1+1)))
}
#dREML
FASTmrEMMA.delta.REML.dLL.c<-function(logdelta,X,M,M.y,yMy,v){
  #X=X_c:n*1,M=M_c:n*n,M.y=M_c%*%y_c:n*1,yMy=t(y_c)%*%M_c%*%y_c:1*1
  #n<-dim(M)[1]
  delta <-  exp(logdelta)
  #ci<-crossprod(X)
  delta1<-as.numeric(t(X)%*%M%*%X)
  xMy<-as.numeric(crossprod(X,M.y))
  return(-0.5*delta1/(1+delta*delta1)+0.5*v/((1+delta*delta1)*(as.numeric(yMy)*(1+delta*delta1)/(xMy^2)-delta)))
}
####################
#20190906
FASTmrEMMA.MLE.c<-function(X,M,M.y,yMy,n, ngrids=100, llim=-10, ulim=10, esp=1e-10){
  
  logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
  m <- length(logdelta)
  #delta <- exp(logdelta)
  dLL<-FASTmrEMMA.delta.ML.dLL.c(logdelta,X,M,M.y,yMy,n)
  optlogdelta <- vector(length=0)
  optLL <- vector(length=0)
  
  if ( dLL[1] < esp ) {
    optlogdelta <- append(optlogdelta, llim)
    optLL <- append(optLL,FASTmrEMMA.delta.ML.LL.c(llim,X,M,M.y,yMy,n))
  }
  if ( dLL[m-1] > 0-esp ) {
    optlogdelta <- append(optlogdelta, ulim)
    #optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    optLL <- append(optLL, FASTmrEMMA.delta.ML.LL.c(ulim,X,M,M.y,yMy,n))
  }
  for( i in 1:(m-1) )
  {
    if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
    {
      #r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
      #r <- uniroot(FASTmrEMMA.delta.ML.dLL.c,lower = logdelta[i],upper = logdelta[i+1],X,M,M.y,yMy,n)
      r <- uniroot(FASTmrEMMA.delta.ML.dLL.c,c(logdelta[i],logdelta[i+1]),X=X,M=M,M.y=M.y,yMy=yMy,n=n)
      
      optlogdelta <- append(optlogdelta, r$root)
      #optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
      optLL <- append(optLL,FASTmrEMMA.delta.ML.LL.c(r$root,X,M,M.y,yMy,n))
    }
  }
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  optLL=replaceNaN(optLL)
  maxLL <- max(optLL)
  xMy<-crossprod(X,M.y)
  xMx<-crossprod(X,(M%*%X))
  maxve <-(yMy-maxdelta*(xMy)^2/(1+maxdelta*xMx))/n
  #(sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
  maxvg <- maxve*maxdelta
  #alpha<-inv()
  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,delta1=xMx,xMy=xMy))
}

FASTmrEMMA.REMLE.c<-function(X,M,M.y,yMy,v, ngrids=100, llim=-10, ulim=10, esp=1e-10){
  
  logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
  m <- length(logdelta)
  #delta <- exp(logdelta)
  dLL<-FASTmrEMMA.delta.REML.dLL.c(logdelta,X,M,M.y,yMy,v)
  optlogdelta <- vector(length=0)
  optLL <- vector(length=0)
  if ( dLL[1] < esp ) {
    optlogdelta <- append(optlogdelta, llim)
    optLL <- append(optLL,FASTmrEMMA.delta.REML.LL.c(llim,X,M,M.y,yMy,v))
  }
  if ( dLL[m-1] > 0-esp ) {
    optlogdelta <- append(optlogdelta, ulim)
    #optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    optLL <- append(optLL, FASTmrEMMA.delta.REML.LL.c(ulim,X,M,M.y,yMy,v))
    
  }
  for( i in 1:(m-1) )
  {
    if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
    {
      #r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
      r <- uniroot(FASTmrEMMA.delta.REML.dLL.c,lower = logdelta[i],upper = logdelta[i+1],X=X,M=M,M.y=M.y,yMy=yMy,v=v)
      
      optlogdelta <- append(optlogdelta, r$root)
      #optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
      optLL <- append(optLL,FASTmrEMMA.delta.REML.LL.c(r$root,X,M,M.y,yMy,v))
    }
  }
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  optLL=replaceNaN(optLL)
  maxLL <- max(optLL)
  xMy<-crossprod(X,M.y)
  xMx<-crossprod(X,(M%*%X))
  maxve <-(yMy-maxdelta*(xMy)^2/(1+maxdelta*xMx))/v
  #(sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
  maxvg <- maxve*maxdelta
  #alpha<-inv()
  
  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,delta1=xMx,xMy=xMy))
}

emma.maineffects.B<-function(Z=NULL,K,deltahat.g,complete=TRUE){
  if( is.null(Z) ){
    return(emma.maineffects.B.Zo(K,deltahat.g))
  }
  else{
    return(emma.maineffects.B.Z(Z,K,deltahat.g,complete))
  }
}

emma.maineffects.B.Zo <-function(K,deltahat.g){
  t <- nrow(K)
  stopifnot(ncol(K) == t)
  B<-deltahat.g*K+diag(1,t)
  eig<-eigen(B,symmetric=TRUE)
  qr.B<-qr(B)
  q<-qr.B$rank
  stopifnot(!is.complex(eig$values))
  A<-diag(1/sqrt(eig$values[1:q]))
  Q<-eig$vectors[,1:q]
  C<-Q%*%A%*%t(Q)
  return(list(mC=C,Q=Q,A=A))
}

emma.maineffects.B.Z <- function(Z,K,deltahat.g,complete=TRUE){
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)  
  B <- deltahat.g*Z%*%K%*%t(Z)+diag(1,n)
  eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
  qr.B<-qr(B)
  q<-qr.B$rank
  stopifnot(!is.complex(eig$values))
  A<-diag(1/sqrt(eig$values[1:q]))
  Q<-eig$vectors[,1:q]
  C<-Q%*%A%*%t(Q)
  return(list(mC=C,Q=Q,A=A,complete=TRUE))
}

emma.MLE0.c <- function(Y_c,W_c){
  n <- length(Y_c)
  stopifnot(nrow(W_c)==n)
  M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
  etas<-crossprod(M_c,Y_c)
  LL <- 0.5*n*(log(n/(2*pi))-1-log(sum(etas*etas)))
  return(list(ML=LL,M=M_c,n=n))
  
}



emma.REMLE0.c <- function(Y_c,W_c){#20190831
  
  n <- length(Y_c)
  stopifnot(nrow(W_c)==n)
  
  t <-qr(W_c)$rank
  v <-n-t
  
  M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
  etas<-crossprod(M_c,Y_c)
  
  LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
  return(list(REML=LL,M=M_c,v=v))
  
}

replaceNaN<-  function(LL) {
  index=(LL=="NaN")
  if(length(index)>0) theMin=min(LL[!index])
  if(length(index)<1) theMin="NaN"
  LL[index]=theMin
  return(LL)    
}


maf.fun<-function(snp){
  leng<-length(snp)
  id.1<-length(which(snp==1))
  id.0<-length(which(snp==0))
  id.0.5<-length(which(snp==0.5))
  maf.1<-id.1/leng
  maf.0.5<-id.0.5/leng
  maf.0<-id.0/leng
  ma1<-(2*id.1+id.0.5)/(2*leng)
  ma2<-(2*id.0+id.0.5)/(2*leng)
  maf.min<-min(ma1,ma2)
  return(list(maf.1,maf.0,maf.0.5,maf.min))
}


pve.fun<-function(beta,maf){
  pve<-(maf$p1-maf$p1^2+0.25*maf$p3-0.25*maf$p3^2-maf$p1*maf$p3)*beta^2
  return(pve)
}

yraw<-matrix(phe[,1],,1)
xnames<-gen[,1:2]
snp1<-gen[,3:ncol(gen)]
mydata<-t(matrix(snp1,nrow=dim(snp1)[1]))
m<-dim(mydata)[2]
n<-dim(mydata)[1]
Y<-yraw
K<-matrix(kk,nrow=dim(kk)[1])
W0<-matrix(1,n,1)

if(is.null(psmatrix)==FALSE){
  W1<-psmatrix
  W<-cbind(W0,W1)
  
}
if(is.null(psmatrix)==TRUE){
  W<-W0 
  
}

rm(kk)
gc()

if(flagREMLE==1){
  remle1<-emma.REMLE(Y, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
}else{
  remle1<-emma.MLE(Y, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
  
}

remle1.deltahat.g<-remle1$delta
remle1.B1<-emma.maineffects.B(Z=NULL,K,remle1.deltahat.g)
C2<-remle1.B1$mC

rm(remle1.B1)
gc()

if(flagREMLE==1){
  
  ys=Y;xs=mydata;Z=C2;X0=W;ngrids=100;llim=-10;ulim=10;esp=1e-10
  ys <- Z%*%ys  
  xs <- Z%*%xs
  X0 <- Z%*%X0
  
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  
  n <- nrow(ys)
  t <- ncol(xs)
  q<- if ( is.matrix(X0) ) ncol(X0) else 1
  v<-n-q
  
  MLE0<-emma.REMLE0.c(ys,X0)
  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  vgs <- vector(length=t)
  ves <- vector(length=t)
  lambdas <- vector(length=t)
  bhats<-vector(length=t)
  d <- vector(length=t)
  
  stats <- vector(length=t)
  ps <- vector(length=t)
  M<-MLE0$M
  M.y<-M%*%ys
  yMy<-crossprod(ys,M.y)
  
  
  cl.cores <- detectCores()
  if ((cl.cores<=2)||(is.null(CLO)==FALSE)){
    cl.cores<-1
  }else if(cl.cores>2){
    if(cl.cores>10){
      cl.cores<-10
    }else {  
      cl.cores <- detectCores()-1
    }
  } 
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl)
  
  REML.LRT.c2<-foreach(i=1:t,.combine = 'rbind')%dopar%{
    
    #MLE1 <- emma.REMLE.c (ys, x0v, K=1, xv, qr.X0,ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)#20181112
    MLE1 <- FASTmrEMMA.REMLE.c(X=xs[,i],M,M.y,yMy,v, ngrids=100, llim=-10, ulim=10, esp=1e-10)
    
    if(is.na(MLE1$REML)==TRUE){
      ps[i]<-1
    }else{
      ML1s[i]<-MLE1$REML
      ML0s[i]<-MLE0$REML
      vgs[i]<-MLE1$vg
      ves[i]<-MLE1$ve
      lambdas[i] <- MLE1$delta
      ###################
      d[i] <- 1/(1+MLE1$delta*MLE1$delta1)
      #bhats[i]<-MLE1$lambda*MLE1$xMy/(1+MLE1$lambda*MLE1$delta1)
      #bhats[i]<-MLE1$delta*MLE1$xMy/d[i]
      bhats[i]<-MLE1$delta*MLE1$xMy*d[i]
      #to record me=sum(d) 
      
      stats[i]<- 2*(MLE1$REML-MLE0$REML)
      ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2
    }
    
    c(ps[i],bhats[i],lambdas[i],d[i],ML1s[i],ML0s[i],stats[i],vgs[i],ves[i])
  }
  stopCluster(cl) 
}else{
  
  #FASTmrEMMA.ML.LRT.c <- function(ys, xs, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
  #20190910
  #Z=C,X0=W=W0,xs=x:snp,n*p
  ys=Y;xs=mydata;Z=C2;X0=W;ngrids=100;llim=-10;ulim=10;esp=1e-10
  ys <- Z%*%ys   
  xs <- Z%*%xs
  X0 <- Z%*%X0
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  n <- nrow(ys)
  t <- ncol(xs)
  q<- if ( is.matrix(X0) ) ncol(X0) else 1
  v<-n-q
  MLE0<-emma.MLE0.c(ys,X0)
  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  vgs <- vector(length=t)
  ves <- vector(length=t)
  lambdas<-vector(length=t)
  bhats<-vector(length=t)
  #   
  d <- vector(length=t)
  stats <- vector(length=t)
  ps <- vector(length=t)
  #n<-199
  #M<-diag(1,n)-X0%*%ginv(crossprod(X0))%*%t(X0)
  M<-MLE0$M
  M.y<-M%*%ys
  yMy<-crossprod(ys,M.y)
  
  cl.cores <- detectCores()
  if((cl.cores<=2)||(is.null(CLO)==FALSE)){
    cl.cores<-1
  }else if(cl.cores>2){
    if(cl.cores>10){
      cl.cores<-10
    }else {  
      cl.cores <- detectCores()-1
    }
  }   
  
  cl <- makeCluster(cl.cores)
  registerDoParallel(cl)
  
  REML.LRT.c2<-foreach(i=1:t,.combine = 'rbind')%dopar%{
    
    MLE1 <- FASTmrEMMA.MLE.c(X=xs[,i],M,M.y,yMy,n, ngrids=100, llim=-10, ulim=10, esp=1e-10)
    if(length(MLE1$vg)!=0){
      ML1s[i]<-MLE1$ML
      ML0s[i]<-MLE0$ML
      vgs[i]<-MLE1$vg
      ves[i]<-MLE1$ve
      lambdas[i]<-MLE1$delta
      ###################
      d[i] <- 1/(1+MLE1$delta*MLE1$delta1)
      #bhats[i]<-MLE1$lambda*MLE1$xMy/(1+MLE1$lambda*MLE1$delta1)
      #bhats[i]<-MLE1$delta*MLE1$xMy/d[i]
      bhats[i]<-MLE1$delta*MLE1$xMy*d[i]
      #to record me=sum(d) 
      stats[i]<- 2*(MLE1$ML-MLE0$ML)
      ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2#20160619
    }else{
      ps[i]<-1
    }             
    c(ps[i],bhats[i],lambdas[i],d[i],ML1s[i],ML0s[i],stats[i],vgs[i],ves[i])
  }
  stopCluster(cl)
} 

rm(Z,xs)
gc()

REML.LRT.c2.new<-data.frame(REML.LRT.c2)

rm(C2,mydata)
gc()

parms<-data.frame(chr.locus=xnames,REML.LRT.c2.new)
names(parms)<-NULL
parms<-as.matrix(parms)


parmeter<-parms[,1:4]
parmeter[,3]<--log10(parmeter[,3])
parmeter[which(abs(parmeter)>1e-4)]<-round(parmeter[which(abs(parmeter)>1e-4)],4)
parmeter[which( abs(parmeter)<1e-4)]<-as.numeric(sprintf("%.4e", parmeter[which( abs(parmeter)<1e-4)]))


if(inputform==1){
  parmsShow<-cbind(genRaw[-1,1],parmeter,genRaw[-1,4])
  parmsShow<-parmsShow[,c(1,2,3,5,4,6)]
  colnames(parmsShow)<-c("Marker","Chromosome","Marker position (bp)","SNP effect (FASTmrEMMA)","'-log10(P) (FASTmrEMMA)'","Genotype for code 1")
  
}
if(inputform==2){
  outATCG<-matrix(outATCG,,1)
  meadd<-matrix(" ",nrow(parms),1)
  parmsShow<-cbind(genRaw[-1,1],parmeter,outATCG)
  parmsShow<-parmsShow[,c(1,2,3,5,4,6)]
  colnames(parmsShow)<-c("Marker","Chromosome","Marker position (bp)","SNP effect (FASTmrEMMA)","'-log10(P) (FASTmrEMMA)'","Genotype for code 1")
  
}
if(inputform==3){
  outATCG<-matrix(outATCG,,1)
  outATCG<-unlist(strsplit(outATCG,""))
  outATCG<-matrix(outATCG[c(TRUE,FALSE)],,1)
  
  parmsShow<-cbind(genRaw[-1,1],parmeter,outATCG)
  parmsShow<-parmsShow[,c(1,2,3,5,4,6)]
  colnames(parmsShow)<-c("Marker","Chromosome","Marker position (bp)","SNP effect (FASTmrEMMA)","'-log10(P) (FASTmrEMMA)'","Genotype for code 1")
}


Xemma<-data.frame(chr.locus=xnames,REML.LRT.c2.new)
vid<-which(as.numeric(Xemma[,3])<=svpal)

if(length(vid)!=0){
  if(length(vid)==1){
    snp.emma.opt<-matrix(gen[vid,],1,)
    xname.emma.opt<-matrix(snp.emma.opt[,1:2],1,)
    #mafall4.opt<-mafall4[vid,]
    snp4<-matrix(snp.emma.opt[,3:dim(snp.emma.opt)[2]],1,)
    xdata<-t(snp4)
    xdata<-matrix(xdata,,1)
  }else{
    snp.emma.opt<-as.matrix(gen[vid,])
    xname.emma.opt<-snp.emma.opt[,1:2]
    #mafall4.opt<-mafall4[vid,]
    snp4<-snp.emma.opt[,3:dim(snp.emma.opt)[2]]
    xdata<-t(snp4)
  }
  
  xdata<-t(snp4)
  ydata<-Y
  u1<-ebayes_EM(x=W,z=xdata,y=ydata)
  emma.lod<-likelihood(xxn=W,xxx=xdata,yn=ydata,bbo=u1$u)
  idslod<-which(emma.lod>=svmlod)
  
  if(length(idslod)!=0){
    # maf.snp.4<-mafall4.opt[idslod,]
    if(length(idslod)==1){
      chrlocus<-matrix(xname.emma.opt[idslod,],1,)
    }else{
      chrlocus<-as.matrix(xname.emma.opt[idslod,])
    }
 
    maf.snp<-apply(snp1,1,maf.fun)
    rm(snp1)
    gc()
    
    maf.snp.1<-unlist(maf.snp)
    maf.snp.2<-matrix(maf.snp.1,nrow = 4)
    maf.snp.3<-t(maf.snp.2)
    maf.snp.4<-data.frame(maf.snp.3)
    names(maf.snp.4)<-c("p1","p2","p3","maf")
  
    pve.opt.all.1<-pve.fun(u1$u[idslod],maf.snp.4[vid,][idslod,])
    pve.opt.all<-pve.opt.all.1/as.vector(max(var(Y),sum(pve.opt.all.1)+u1$sigma2))*100
    
    qtneffect<-matrix(u1$u[idslod],,1)
    lodscore<-matrix(emma.lod[idslod],,1)
    log10P <- as.matrix(-log10(pchisq(lodscore*4.605,1,lower.tail = F)))
    maff<-matrix(maf.snp.4[vid,][idslod,]$maf,,1)
    r2<-matrix(pve.opt.all,,1)
    wanbefore<-cbind(qtneffect,lodscore,log10P,r2,maff)
    
    wanbefore[which(abs(wanbefore)>1e-4)]<-round(wanbefore[which(abs(wanbefore)>1e-4)],4)
    wanbefore[which(abs(wanbefore)<1e-4)]<-as.numeric(sprintf("%.4e", wanbefore[which(abs(wanbefore)<1e-4)]))
    wanbefore <- matrix(wanbefore,,5)
    
    wan<-cbind(chrlocus,wanbefore)
    phenotype.var<-var(Y)
    sigma2<-u1$sigma2
    pee<-matrix("",dim(wan)[1],1)
    vess<-matrix("",dim(wan)[1],1)
    pee[1]<-round(phenotype.var,4)
    vess[1]<-round(sigma2,4)
    
    if(inputform==1){
      genRaw<-as.data.frame(genRaw)
      genraw<-genRaw[-1,1:4]
      
      wan_len<-dim(wan)[1]
      marker<-character()
      snp<-character()
      
      for(i in 1:wan_len){
        chr_pos<-which(genraw[,2]==wan[i,1])
        new_matrix<-genraw[chr_pos,]
        posi_pos<-which(new_matrix[,3]==wan[i,2])
        mark<-matrix(new_matrix[posi_pos,1],1,)
        marker<-rbind(marker,mark)
        sn<-matrix(new_matrix[posi_pos,4],1,)
        snp<-rbind(snp,sn)
      }
    }
    if(inputform==2){
      genRaw<-as.data.frame(genRaw)
      genraw<-genRaw[-1,1:4]
      wan_len<-dim(wan)[1]
      marker<-character()
      snp<-character()
      for(i in 1:wan_len){
        chr_pos<-which(genraw[,2]==wan[i,1])
        new_matrix<-genraw[chr_pos,]
        posi_pos<-which(new_matrix[,3]==wan[i,2])
        mark<-matrix(new_matrix[posi_pos,1],1,)
        marker<-rbind(marker,mark)
        sn<-matrix(new_matrix[posi_pos,4],1,)
        snp<-rbind(snp,sn)
      }
      
    }
    if(inputform==3){
      genRaw<-as.data.frame(genRaw)
      genraw<-genRaw[-1,c(1,3,4,12)]
      
      wan_len<-dim(wan)[1]
      marker<-character()
      snp<-character()
      for(i in 1:wan_len){
        chr_pos<-which(genraw[,2]==wan[i,1])
        new_matrix<-genraw[chr_pos,]
        posi_pos<-which(new_matrix[,3]==wan[i,2])
        mark<-matrix(new_matrix[posi_pos,1],1,)
        marker<-rbind(marker,mark)
        sn<-matrix(new_matrix[posi_pos,4],1,)
        snp<-rbind(snp,sn)
      }
    }
    
    wan<-cbind(marker,wan,snp,vess,pee)
    colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen(total)")
    wan<-as.data.frame(wan)
    
  }  
}  
output<-list(result1=parmsShow,result2=wan)
return(output)
}
}

