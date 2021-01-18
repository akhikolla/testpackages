###############################################################################
##
##   Author : Veronique CARIOU (ONIRIS) veronique.cariou@oniris-nantes.fr
##
#############################################################library(plyr)
nmodeMatrix<-function(X,Mode=1) {
  dimX=dim(X);
  ncols<-prod(dimX[-Mode])
  permode<-matrix(0,nrow=3,ncol=3)
  permode[1,]=c(1, 2 ,3);
  permode[2,] =c(2, 1 ,3);
  permode[3,]=c(3,1, 2 );
  perX=aperm(X,permode[Mode,]);
  mat<-matrix(as.vector(perX),nrow=dimX[Mode],ncol=ncols)
  return(mat)
}

Product3D1D<-function(X,vec,Mode) {
  v = dim(X);
  M<-nmodeMatrix(X,Mode=Mode);
  XX <- matrix(as.vector(vec%*%M),nrow=v[-Mode][1],ncol=v[-Mode][2]);
  return(XX)
}




##Estimation of the one-rank Candecomp / Parafac model
#start=0 : initial solution given by SVD
#start=1 : inital solution is random
#start=2 : initial solution (mode1 & mode 2) is given
#start=3 : intial solution given by SVD of the mean
#NN      : Non negativity constraint TRUE or FALSE
CP1<-function(X,NN=FALSE, start=1,initMode1=NULL,initMode2=NULL,initMode3=NULL,itermax=100000,tol=1e-06){
  if (start==0) {
    M1<-nmodeMatrix(X,1)
    u1<-svd(M1)$u[,1]
    M2<-nmodeMatrix(X,2)
    v1<-svd(M2)$u[,1]
    if (NN) {
      if (length(which(v1<0))) {
        if (length(which(v1<0))==length(v1)) {
          u1<--u1
          v1<--v1
        }
        else {
          v1[which(v1<0)]<-0
        }
      }
    }
  }
  else {
    if (start==1) {
      u1<-runif(dim(X)[1],0,1)
      u1<-u1/norm(as.matrix(u1))
      v1<-runif(dim(X)[2],0,1)
      v1<-v1/norm(as.matrix(v1))
    }
    else {
      if (start==2) {
        u1<-initMode1
        u1<-u1/norm(as.matrix(u1))
        v1<-initMode2
        v1<-v1/norm(as.matrix(v1))
      }
      else {
        M<-plyr::aaply(.data=X,.margins=c(1,2),.fun=mean)
        res<-svd(M)
        u1<-res$u[,1]
        v1<-res$v[,1]
        if (NN) {
          if (length(which(v1<0))) {
            if (length(which(v1<0))==length(v1)) {
              u1<--u1
              v1<--v1
            }
            else {
              v1[which(v1<0)]<-0
            }
          }
        }
      }
    }
  }
  iter<-1
  deltafit<-1
  fit<-NULL
  normX<-sum(X^2)
  fiti<-sum(X^2)
  while ((iter <itermax) & (deltafit>tol)) {
    M3<-nmodeMatrix(X,3)
    w1<-(M3%*%kronecker(v1,u1,FUN="*"))/(sum(v1^2)*sum(u1^2))
    M2<-nmodeMatrix(X,2)
    v1<-(M2%*%kronecker(w1,u1,FUN="*"))/(sum(w1^2)*sum(u1^2))
    ineg<-which(v1<0)
    ipos<-which(v1>=0)
    if (length(ipos)==0)
    {
      v1 <- -v1
      w1 <- -w1
    }
    else {
      if (length(ineg) && (sum(v1[ineg]^2)>sum(v1[ipos]^2))) {
        v1 <- -v1
        w1 <- -w1
        
      }
    }
    if (NN) {
      ineg<-which(v1<0)
      if (length(ineg)) {
        v1[ineg] <- 0
      }
    }
    M1<-nmodeMatrix(X,1)
    u1<-(M1%*%kronecker(w1,v1,FUN="*"))/(sum(w1^2)*sum(v1^2))
    R1<-array(0,dim=dim(X))
    for( k in 1:dim(X)[3]) {
      R1[,,k]=w1[k]*(u1%o%v1);
    }
    sum.R1.2 <- sum(R1^2)
    deltafit<-abs(fiti-(sum.R1.2/normX))
    #bascule en fit et loss
    fiti<-sum.R1.2/normX;
    fit<-c(fit,fiti)
    iter<-iter+1
  }
  #normalisation des vecteurs
  lambda<-sqrt( sum(u1^2))+ sqrt( sum(v1^2))+ sqrt( sum(w1^2))
  vp<-sum.R1.2
  u1<-u1/sqrt( sum(u1^2))
  v1<-v1/sqrt( sum(v1^2))
  w1<-w1/sqrt( sum(w1^2))
  res<-list(u=u1,v=v1,w=w1,lambda=lambda,loss=(1-fit),fit=fit,afit=fiti*normX,aloss=sum((X-R1)^2),vp=vp)
}


##Estimation of the one-rank Candecomp / parafac model with multiple starts
##given Als procedure
## NN = FALSE : Non Negativity constraint
CP1_MS<-function(X,NN=FALSE,nrandom=10,itermax=1000,tol=1e-06){
  globalres <- NULL
  allsol <- vector("numeric",length=nrandom)
  for (i in 1:nrandom) {
    res<-CP1(X=X,NN=NN,start=1,itermax=itermax,tol=tol)
    allsol[i] <- res$afit
    if ((i==1) || (res$afit>globalres$afit)) {
        globalres<-res
    }
  }
  globalres$allsol <- allsol
  return(globalres)
}

cp.deflate <- function(X,cp1) {
  projector<-cp1$u%*%(t(cp1$u)%*%cp1$u)^(-1) %*%t(cp1$u)
       for (k in 1:dim(X)[[3]]) {
          X[,,k] <- X[,,k] - projector %*% X[,,k];
       }
  return(X)
}
## Estimation of a two-rank candecomp parafac model
## given a greedy parafac algorithm

CP2_MS<-function(X,ncp=2,nrandom=10,itermax=1000,tol=1e-06){
  U <- matrix(0,nrow=dim(X)[1],ncol=ncp)
  V <- matrix(0,nrow=dim(X)[2],ncol=ncp)
  W <- matrix(0,nrow=dim(X)[3],ncol=ncp)
  lambda <- vector("numeric",length=ncp)
  afit <-vector("numeric",length=ncp)
  aloss <- vector("numeric",length=ncp)
  for (k in 1:ncp) {
    res.CP <- CP1_MS(X=X,NN=FALSE,nrandom=nrandom,itermax=itermax,tol=tol)
    U[,k] <- res.CP$u
    V[,k] <- res.CP$v
    W[,k] <- res.CP$w
    lambda[k] <- res.CP$lambda
    afit[k] <- res.CP$afit
    aloss[k] <- res.CP$aloss
    #deflation
    X <- cp.deflate(X,res.CP)
  }
  globalres <- list(u=U,v=V,w=W,lambda=lambda,afit=afit,aloss=aloss)
  return(globalres)
}
