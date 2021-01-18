MINLDET <- -500

CheckMSing <- function(M,limlnk2,limdiag=NULL,scale,fullcheck=TRUE)
{
   if ( any(!is.finite(M)) || (!is.null(limdiag) && any(diag(M)<limdiag)) ) return(TRUE)
   if (scale) M <- cov2cor(M)
   trM <- sum(diag(M)) 
   lnMdet <- determinant(M,logarithm=TRUE)
   if (lnMdet$sign==-1) return(TRUE) 
   p <- ncol(M)
   if ( p*log(trM) - lnMdet$modulus < limlnk2) return(FALSE)   
   if (fullcheck) {   
     Megval <- eigen(M,symmetric=TRUE,only.values=TRUE)$values
     if (Megval[p] <= 0.) return(TRUE)
     if (log(Megval[1]/Megval[p]) > limlnk2) return(TRUE)
   }
   FALSE  # return(FALSE)
}

CheckSigmaSing <- function(Cf,Sigma,limlnk2,scale,limdiag=1e-300)
{
    if (any(diag(Sigma)<limdiag)) return(TRUE)
    if (scale) Sigma <- cov2cor(Sigma)
    if (Cf==1) return(CheckMSing(Sigma,limlnk2,scale=FALSE))
    p <- ncol(Sigma)
    q <- p/2
    maxegval <- 0.
    minegval <- Inf

    if (Cf==2) {
      for (j in 1:q) {
        a <- Sigma[j,j]
        b <- Sigma[q+j,q+j]
        c <- Sigma[j,q+j]
        trSigj <- a + b
        detSigj <- a*b - c^2
        d <- sqrt(trSigj^2-4*detSigj)
        maxegval <- max(maxegval,trSigj+d)
        minegval <- min(minegval,trSigj-d)
        if (minegval<=0.) return(TRUE) 
      }
      if (log(maxegval/minegval) > limlnk2) return(TRUE)
      return(FALSE)
    }

    if (Cf==3) {
      if (CheckMSing(Sigma[1:q,1:q],limlnk2,scale=FALSE,fullcheck=FALSE)) return(TRUE)
      if (CheckMSing(Sigma[(q+1):p,(q+1):p],limlnk2,scale=FALSE,fullcheck=FALSE)) return(TRUE)
      MidPEigv <- eigen(Sigma[1:q,1:q],symmetric=TRUE,only.values=TRUE)$values
      LogREigv <- eigen(Sigma[(q+1):p,(q+1):p],symmetric=TRUE,only.values=TRUE)$values
      maxegval <- max(MidPEigv[1],LogREigv[1])
      minegval <- min(MidPEigv[q],LogREigv[q])
      if (minegval<=0.) return(TRUE) 
      if (log(maxegval/minegval) > limlnk2) return(TRUE)
      return(FALSE)
    }

    if (Cf==4) {
      vars <- diag(Sigma) 
      egval1 <- max(vars)
      egvalp <- min(vars)
      if (egvalp < 0 || log(egval1/egvalp) >=  limlnk2) return(TRUE)
      else return(FALSE)
    }
}

CheckSigmak <- function(Cf,Sigmak,MaxVarGRt,limlnk2,scale)
{
   p <- dim(Sigmak)[1]
   k <- dim(Sigmak)[3]
   minVar <- rep(Inf,p)
   maxVar <- rep(0.,p)
   for (g in 1:k) { 
     if (CheckSigmaSing(Cf,Sigmak[,,g],limlnk2,scale=scale)==TRUE) return(list("Singular"))
     gVar <- diag(Sigmak[,,g])
     minVar <- pmin(minVar,gVar)
     maxVar <- pmax(maxVar,gVar)
   }
   VarGRt <- maxVar/minVar
   viollst <- which(VarGRt > MaxVarGRt)
   nviol <- length(viollst) 
   if (nviol > 0) {
     mingrpVar <- integer(nviol)
     for (v in 1:nviol) mingrpVar[v] <- which.min(Sigmak[v,v,])    
     return(list("LargeVarR",viollst=viollst,mingrpVar=mingrpVar,Maxviol=max(VarGRt[viollst])/MaxVarGRt))
   }
   return(list("Valid"))
}

safecholR <- function(Sigma)
{
   R <- cov2cor(Sigma)
   if ( any(!is.finite(R)) || any(abs(R)>1.) ) return(diag(1.,nrow(Sigma)))
   chol(R)
}

Safepdsolve <- function(M, maxlnk2, scale, limdiag=1e-100)
{
  if (scale) {  
    p <- ncol(M)
#    diagsrqt <- numeric(p)
    Mscl <- matrix(nrow=p,ncol=p)
#    for (i in 1:p) diagsrqt[i] <- sqrt(M[i,i])
    dM <- diag(M)
    if (any(dM < 0.)) return(NULL)
    diagsrqt <- sqrt(dM)
    for (i in 1:p) for (j in 1:i) Mscl[i,j] <- Mscl[j,i] <- diagsrqt[i]*diagsrqt[j]
    M <- M / Mscl
  }
  if ( CheckMSing(M,limlnk2=maxlnk2,limdiag=limdiag,scale=FALSE) )  return(NULL)
  MInv <- solve(M)
  if (scale) return(MInv / Mscl)
  else return(MInv) 
}

VectGaussLogLik <- function(x,p,u,SigmaInv,CheckSing,k2max,cnst=NULL)
{

#  Computes the log-likelihood of the mean, u and covariance matrix, Sigma estimates, 
#  for a multivariate normal distribution, given an  observed vector, x 
#  When CheckSing is set to TRUE it returns a "large penalty" if Sigma is found to be numerically singular
#
#  Parameters:
#  
#  p          -  x dimension
#  x          -  vector of observation data
#  u          -  vector of mean estimates 
#  SigmaInv   -  The inverse of the covariance matrix estimate 
#  k2max      -  maximum (L2 norm) condition number for Sigma to be considered numerically non-singular  
#  CheckSing  -  Boolean flag that should be set to TRUE when it is necessary to check the singularity of Sigma
#  cnst       -  the log-likeklihood constant -1/2 ln(2 Pi |Sigma|), when already known (or null if needs to be computed)

#   value - the log-likelihood of u and Sigma, given x

  if (CheckSing==TRUE)  {
    PenF <- 1E6
    egvalues <- eigen(SigmaInv,symmetric=TRUE,only.values=TRUE)$values
    k2 <- egvalues[p]/abs(egvalues[1])
    if (!is.finite(k2)) return(-PenF*k2max)
    if (k2 > k2max) return(PenF*(k2max-k2))
  }

  if (is.null(cnst))  {
    lndetSig <- 1./determinant(SigmaInv,logarithm=TRUE)$modulus
    cnst <- -(p*log(2*pi)+lndetSig)/2
  }
  dev <- x-u

  cnst -(dev%*%SigmaInv%*%dev)/2
}  

MDataGaussLogLik <- function(n,p,X,u,Sigma=NULL,SigmaInv=NULL,cnst=NULL,k2max,invalcode=NULL)
{

#  Computes the log-likelihood of the mean, u and covariance matrix, Sigma estimates, 
#  for a multivariate normal distribution, given an  observed data matrix, X 
#  When Sigma is found to be numerically singular returns a "large penalty" 
#  if invalcode is set to NULL, or invalcode otherwise 
#
#  Parameters:
#  
#  p          -  x dimension
#  X          -  matrix of observation data
#  u          -  vector of mean estimates 
#  Sigma      -  The covariance matrix estimate 
#  k2max      -  maximum (L2 norm) condition number for Sigma to be considered numerically non-singular  

#   value - a vector with  the log-likelihood contribution of each row of the X data matrix

  minlndet <- -500    # minimum log-determinant for Sigma to be considered numerically non-singular 
  PenF <- 1e6	      # penalty factor for numerically singular covariance matrices	

  if (is.null(cnst)) {
    lndetSig <- determinant(Sigma,logarithm=TRUE)$modulus
    if (!is.finite(lndetSig) || lndetSig < minlndet) {
      if (is.null(invalcode)) return(rep(PenF*(lndetSig-minlndet),n))
      else return(rep(invalcode,n))
    }
    cnst <- -(p*log(2*pi)+lndetSig)/2
  }
  if (is.null(SigmaInv)) {
    SigmaInv <- Safepdsolve(Sigma,maxlnk2=log(k2max),scale=TRUE)
    if (is.null(SigmaInv)) {
      if (is.null(invalcode)) return(rep(PenF,n))
      else return(rep(invalcode,n))
    }
  } 

  apply(X,1,VectGaussLogLik,p=p,u=u,SigmaInv=SigmaInv,CheckSing=FALSE,k2max=NULL,cnst=cnst)
}  

MClusLikz <- function(z,X,n,p,k,Homoc,k2max,penalty=-1.e99)
{
  Vnames <- names(X)
  Onames <- rownames(X)
  Likk <- matrix(nrow=n,ncol=k)
  LikExp <- matrix(nrow=n,ncol=k)
  muk <- matrix(nrow=p,ncol=k,dimnames=list(Vnames,NULL))	
  Repmuk <- array(dim=c(n,p,k),dimnames=list(Onames,Vnames,NULL))
  Wk <- array(dim=c(p,p,k))	
  if (Homoc==FALSE) Sigmak <- array(dim=c(p,p,k))	

  for (g in 1:k)  {
    nk <- apply(z,2,sum)
    tau <- nk/n
    muk[,g] <- apply(matrix(rep(z[,g],p),n,p,byrow=FALSE)*X,2,sum)/nk[g]
    Repmuk[,,g] <- matrix(rep(muk[,g],n),n,p,byrow=TRUE)
    wdev <- as.matrix(matrix(rep(sqrt(z[,g]),p),n,p) * (X - Repmuk[,,g]))
    Wk[,,g] <- t(wdev)%*%wdev
    if (Homoc==FALSE) {
      Sigmak[,,g] <- Wk[,,g]/nk[g]
      LikExp[,g] <- MDataGaussLogLik(n,p,X,muk[,g],Sigma=Sigmak[,,g],k2max=k2max,invalcode=penalty) 
    }
  }	
  if (Homoc==TRUE) {
    Sigma <- apply(Wk,c(1,2),sum)/n
    dimnames(Sigma) <- list(Vnames,Vnames) 
    lndetSig <- determinant(Sigma,logarithm=TRUE)$modulus
    if (!is.finite(lndetSig) || lndetSig < MINLDET) return(penalty)
    cnst <- -(p*log(2*pi)+lndetSig)/2
    SigmaInv <- Safepdsolve(Sigma,maxlnk2=log(k2max),scale=TRUE)
    if (is.null(SigmaInv)) return(penalty)
    for (g in 1:k) LikExp[,g] <- MDataGaussLogLik(n,p,X,muk[,g],SigmaInv=SigmaInv,cnst=cnst,k2max=k2max,invalcode=penalty) 
  }
   if (any(!is.finite(LikExp))) return(penalty)
   nrmfct <- apply(LikExp,1,max)
   LikExp <- sweep(LikExp,1,nrmfct)
   for (g in 1:k) Likk[,g] <- tau[g] * exp(LikExp[,g])
   lnLikall <- log(apply(Likk,1,sum)) 
   if (any(!is.finite(lnLikall))) return(penalty)
   sum(nrmfct+lnLikall)  # return(sum(nrmfct+lnLikall))   
}

MClusLikpar <- function(X,n,p,k,tau,muk,Sigma=NULL,Sigmak=NULL,Homoc,k2max,penalty=-1.e99)
{
  Likk <- matrix(nrow=n,ncol=k)
  LikExp <- matrix(nrow=n,ncol=k)

  if (!is.null(Sigma)) {
    lndetSig <- determinant(Sigma,logarithm=TRUE)$modulus
    if (!is.finite(lndetSig) || lndetSig < MINLDET) return(penalty)
    cnst <- -(p*log(2*pi)+lndetSig)/2
    SigmaInv <- Safepdsolve(Sigma,maxlnk2=log(k2max),scale=TRUE)
    if (is.null(SigmaInv)) return(penalty)
    for (g in 1:k) LikExp[,g] <- MDataGaussLogLik(n,p,X,muk[,g],SigmaInv=SigmaInv,cnst=cnst,k2max=k2max,invalcode=penalty) 
  } else if (!is.null(Sigmak)) {
    for (g in 1:k)  LikExp[,g] <- MDataGaussLogLik(n,p,X,muk[,g],Sigma=Sigmak[,,g],k2max=k2max,invalcode=penalty) 
  } else {
    stop("Arguments Sigma and Sigmak cannot be simultaneusly NULL\n")
  }
   if (any(!is.finite(LikExp))) return(penalty)
   nrmfct <- apply(LikExp,1,max)
   LikExp <- sweep(LikExp,1,nrmfct)
   for (g in 1:k) Likk[,g] <- tau[g] * exp(LikExp[,g])
   lnLikall <- log(apply(Likk,1,sum)) 
   if (any(!is.finite(lnLikall))) return(penalty)
   sum(nrmfct+lnLikall)  # return(sum(nrmfct+lnLikall))   
}

#Addgrp <- function(X,prvz,Cf,Homoc,k2max,tautol,nnewgrps=1)
Addgrp <- function(X,prvz,Cf,Homoc,k2max,tautol,MaxVarGRt,nnewgrps=1)
{
  n <- nrow(X)
  p <- ncol(X)
  prvk <- ncol(prvz)
  newk <- prvk+nnewgrps
  newz <- matrix(nrow=n,ncol=newk,dimnames=list(rownames(prvz),c(colnames(prvz),paste("CP",(prvk+1):newk,sep=""))))
  newz[,1:prvk] <- prvz*(1-nnewgrps*tautol)
  newz[,(prvk+1):newk] <- tautol
#  newsol <- MStep(X,newz,Cf,Homoc,k2max,tautol)
  newsol <- MStep(X,newz,Cf,Homoc,k2max,tautol,MaxVarGRt)
  if (!is.null(newsol)) {
    newsol$z <- newz
    if (Homoc==TRUE) newsol$LnLik <- MClusLikpar(X,n,p,newk,tau=newsol$tau,muk=newsol$muk,Sigma=newsol$Sigma,Sigmak=NULL,Homoc=TRUE,k2max=k2max)
    else newsol$LnLik <- MClusLikpar(X,n,p,newk,tau=newsol$tau,muk=newsol$muk,Sigma=NULL,Sigmak=newsol$Sigmak,Homoc=FALSE,k2max=k2max)
  }

  newsol # return(newsol)
}

