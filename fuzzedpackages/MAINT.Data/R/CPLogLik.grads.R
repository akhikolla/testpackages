library(sn)

#Note:  All these functions assume that vectors with non-diferent elements of symmetrical
#       matrices use a lower-triangular format

ltdind <- function(p) {
  if (p==1) return(1)
  i <- (p-1):1 
  c(1,(p*(p+1)-i*(i+1))/2+1)
}
utdind <- function(p) {i <- 1:p ;  i*(i+1)/2}

fmdind <- function(p) (1:p)^2

dvecind <- function(p)
{
  ind <- 1:p
  if (p>1) for (r in 2:p) ind <- c(ind,(r-1)*p+(r:p))
  ind
}

dvec0ind <- function(p)
{
  ind <- 2:p
  if (p>2) for (r in 2:(p-1)) ind <- c(ind,(r-1)*p+((r+1):p))
  ind
}

lftD <- function(vM,p) 
{
  tmpM <- matrix(nrow=p,ncol=p)
  tmpM[row(tmpM)>=col(tmpM)] <- vM 
  tmpM[row(tmpM)<col(tmpM)] <- t(tmpM)[row(tmpM)<col(tmpM)]
  matrix(tmpM,p^2,1)
}  

lftDplus <- function(vecM,p,indices=dvecind(p),nrows=p*(p+1)/2) 
  matrix(vecM[indices],nrows,1)  

rghtD <- function(vecM,p,ncols=p*(p+1)/2) 
{
  tmpM <- matrix(vecM,p,p)
  tmpM[row(tmpM)>col(tmpM)] <- 
  tmpM[row(tmpM)>col(tmpM)] + t(tmpM)[row(tmpM)>col(tmpM)]  
  matrix(tmpM[row(tmpM)>=col(tmpM)],1,ncols)
}         

#OmgInv.grad <- function(x,p)
OmgInv.grad <- function(x,p,limlnk2)
{
  Omega <- matrix(nrow=p,ncol=p)
  Omega[lower.tri(Omega,diag=TRUE)] <- x
  Omega[upper.tri(Omega)] <- t(Omega)[upper.tri(Omega)] 
#  OmegaI <- pdwt.solve(Omega)
#  OmegaI <- pdwt.solve(Omega, silent=TRUE)
  OmegaI <- Safepdsolve(Omega,maxlnk2=limlnk2,scale=TRUE)
  if(is.null(OmegaI)) {
    warning("Singular Omega matix found in the course of the mle estimation\n")
    return(matrix(0.,nrow=nvcovpar,ncol=nvcovpar)) 
  }  
  
  nvcovpar <- p*(p+1)/2
  Jacob <- matrix(nrow=nvcovpar,ncol=nvcovpar)
  ind <- 0
  for (r in 1:p) for (c in r:p) {
    ind <- ind+1 
    tmpM <- -outer(OmegaI[,r],OmegaI[,c])
    if (r!=c) tmpM <- tmpM -outer(OmegaI[,c],OmegaI[,r])
    Jacob[,ind] <- tmpM[lower.tri(tmpM,diag=TRUE)]
  }   
  Jacob  
}

#Theta.ll.grad <-function(ksi,Omega,eta,y,n,p,OmegaInv)
Theta.ll.grad <-function(ksi,Omega,eta,y,n,p,OmegaInv,limlnk2)
{
  logDet <- attr(OmegaInv,"log.det")
  if (is.matrix(y)) {
    y0 <- scale(y,center=ksi,scale=FALSE)
    nOmgminusS0 <- n*Omega - t(y0)%*%y0
  }  
  else {
    y0 <- y-ksi
    nOmgminusS0 <- Omega - outer(y0,y0)
  }
  diag(nOmgminusS0) <- diag(nOmgminusS0)/2 
  OmegaInvgrad <- nOmgminusS0[lower.tri(nOmgminusS0,diag=TRUE)] 
  Omegapar <- Omega[lower.tri(Omega,diag=TRUE)]
  Omegagrad <-  drop( OmegaInvgrad %*% OmgInv.grad(Omegapar,p,limlnk2) )
  if (is.matrix(y0)) {
    ksigrad <- drop( OmegaInv%*%matrix(apply(y0,2,sum),p,1) - 
                       sum(sapply(y0%*%eta,zeta1))*eta )
    etagrad <- drop( apply(sapply(y0%*%eta,zeta1)*y0,2,sum) )
  }  
  else  {
    ksigrad <- drop( OmegaInv%*%y0 - zeta(1,y0%*%eta)*eta )
    etagrad <- zeta(1,y0%*%eta)*y0
  }  
  c(ksigrad,Omegagrad,etagrad)
}

msnDP.ll.grad <-function(param,y,n=ifelse(is.matrix(y),nrow(y),1),
                    p=ifelse(is.matrix(y),ncol(y),length(y)),limlnk2)
{
  nvcovpar <- p*(p+1)/2
  ksi <- param[1:p]
  Omega <- matrix(nrow=p,ncol=p)
  Omega[lower.tri(Omega,diag=TRUE)] <- param[(p+1):(p+nvcovpar)]
  Omega[upper.tri(Omega)] <- t(Omega)[upper.tri(Omega)] 
  omega <- sqrt(diag(Omega))
  alpha <- param[p+nvcovpar+1:p]
  eta <- alpha/omega
#  Thetagrad <- Theta.ll.grad(ksi,Omega,eta,y,n,p,pdwt.solve(Omega,log.det=TRUE))
  OmegaI <- Safepdsolve(Omega,maxlnk2=limlnk2,scale=TRUE)
  if(is.null(OmegaI)) {
    warning("Singular Omega matix found in the course of the mle estimation\n")
    return(rep(0.,2*p+nvcovpar)) 
  } else {  
    attr(OmegaI,"log.det") <- determinant(Omega,logarithm=TRUE)$modulus
  }
  Thetagrad <- Theta.ll.grad(ksi,Omega,eta,y,n,p,OmegaI,limlnk2=limlnk2)
  etagrad <- Thetagrad[(p+nvcovpar+1):(2*p+nvcovpar)]
  alphagrad <- etagrad/omega
  Omegagrad <- Thetagrad[(p+1):(p+nvcovpar)]
  diagOmegaind <- ltdind(p)
  Omegagrad[diagOmegaind] <- Thetagrad[p+diagOmegaind] - etagrad*(alpha/omega^3)/2
  
  c(Thetagrad[1:p],Omegagrad,alphagrad)
} 

b <- sqrt(2./pi)
b0 <- 2/(4-pi)
zeta1 <- function(x) zeta(1,x)

msnCP.ll.grad <-function(param,y,n=ifelse(is.matrix(y),nrow(y),1),
  p=ifelse(is.matrix(y),ncol(y),length(y)),inoptm=TRUE,
#  PenF=1E12,ldRtol=log(1E-5)+1-p,c2tol=1E-6)
  PenF=1E12,ldRtol=log(1E-5)+1-p,c2tol=1E-6,limlnk2)
{
  if ( (is.matrix(y) && p!=ncol(y)) || (!is.matrix(y) && p!=length(y)) ) 
  {
    stop("Dimension of y is not compatible with argument p\n")
  }
  if (!all(is.finite(param))) { return(rep(0.,length(param))) }
  nvcovpar <- p*(p+1)/2
    
  mu <- param[1:p]
  Sigma <- matrix(nrow=p,ncol=p)
  Sigma[lower.tri(Sigma,diag=TRUE)] <- param[(p+1):(p+nvcovpar)]
  Sigma[upper.tri(Sigma)]  <- t(Sigma)[upper.tri(Sigma)]
  gamma1 <- param[(p+nvcovpar+1):(2*p+nvcovpar)]
  intres <- CPtoDP.ll.grad(mu,Sigma,gamma1,p,inoptm,PenF,ldRtol,c2tol)
  DPll.grad(intres$ksi,intres$Omega,intres$eta,intres$OmegaInv,
    intres$nvcovpar,intres$penaltygrad,
    intres$D23,intres$D32,intres$D33,intres$Dtld32,intres$Dtld33,
    y,n,p,PenF,ldRtol,c2tol,limlnk2=limlnk2)
}

CPtoDP.ll.grad <-function(mu,Sigma,gamma1,p,inoptm=FALSE,PenF=0.,
# ldRtol=log(1e-5)+1-p,c2tol=1e-6,beta0tol=1e-6,limlnk2)
 ldRtol=-500,c2tol=1e-6,beta0tol=1e-6,limlnk2)
{
  nvcovpar <- p*(p+1)/2
  sigma <- sqrt(diag(Sigma))
  mu0 <- sign(gamma1)*sigma*(b0*abs(gamma1))^(1/3)
  SigmaI <- Safepdsolve(Sigma,maxlnk2=limlnk2,scale=TRUE)
  if (is.null(SigmaI))  {
    if (inoptm)  
    {
      warning("Non-positive definite covariance matrix found during gradient computations (which returned 0).\n")
      return(rep(0.,2*p+nvcovpar))
    }  else {
      return(NULL)
    }  
  } else {  
    attr(SigmaI,"log.det") <- determinant(Sigma,logarithm=TRUE)$modulus
  }
  lRdet <- attr(SigmaI,"log.det") - sum(log(diag(Sigma))) 
  if (lRdet < ldRtol)
  {
    mupgrad <- gamma1pgrad <- rep(0.,p)
    tmpM <- 2*SigmaI
    diag(tmpM) <- diag(tmpM)/2 - 1./diag(Sigma)
    Sigmapgrad <- PenF*(ldRtol-lRdet)*tmpM[row(tmpM)>=col(tmpM)]
    penaltygrad <- c(mupgrad,Sigmapgrad,gamma1pgrad)   
  } else {
    penaltygrad <- rep(0.,2*p+nvcovpar)
  } 
  
  SigmaImu0 <- drop(SigmaI%*%mu0)
  beta02 <- drop(mu0%*%SigmaImu0)
  beta0 <- sqrt(beta02)
  p2 <- p^2

  Ip <- diag(p)
  D23 <- apply(kronecker(Ip,mu0)+kronecker(mu0,Ip),2,
               lftDplus,p=p,nrows=nvcovpar)

  if (beta0 < beta0tol)  {
    Dtld32 <- matrix(0.,p,p*(p+1)/2)
    Dtld33 <- matrix(0.,p,p)
  } else {
    mu0bar <- mu0/(sigma*beta0)
    tmpsum <- matrix(0.,p,p*(p+1)/2)
    for (i in 1:p)  {
      tmp <- rep(0.,p2)
      tmp[(i-1)*p+i] <- mu0bar[i]
      tmpsum[i,] <- tmpsum[i,] + rghtD(tmp,p=p)/sigma[i]
    }
    Dtld32 <- beta0*tmpsum/2
    Dtld33 <- (b0/(3*beta02))*sigma*diag(1./mu0bar^2)
  }

  ksi <- mu - mu0 
  Omega <- Sigma + outer(mu0,mu0) 
#  OmegaInv <- pdwt.solve(Omega,silent=TRUE,log.det=TRUE)
  OmegaInv <- Safepdsolve(Omega,maxlnk2=limlnk2,scale=TRUE)
  if (is.null(OmegaInv)) {
    if (inoptm)  
    {
      warning("Non-positive definite matrix Omega found during gradient computations (which returned 0).\n")
      return(rep(0.,2*p+nvcovpar))
    }  else {
      return(NULL)
    }  
  } else {  
    attr(OmegaInv,"log.det") <- determinant(Omega,logarithm=TRUE)$modulus
  }
  b2 <- b^2
  omega <- sqrt(diag(Omega))
  delta <- mu0/(b*omega)
  sclMat <- 1./outer(omega,omega)
  OmgbI <- OmegaInv / sclMat 
  OmgbIdelta <- drop(OmgbI%*%delta)
  c2 <- 1. - delta %*% OmgbIdelta

  if (c2 < c2tol) {
    
    mugrad <- rep(0.,p)

    dOmgbdOmg <- OmgtoOmgbar.grad(Omega,p,sclMat=sclMat)
    tmpM <- 2*outer(OmgbIdelta,OmgbIdelta)
    diag(tmpM) <- diag(tmpM)/2
    vtmpM <- tmpM[lower.tri(tmpM,diag=TRUE)] 
    vdOmgb <- vtmpM %*% dOmgbdOmg  
    Sigmagrad <- vdOmgb + vdOmgb%*%D23%*%Dtld32 
    
    dAomgm1 <- (1/3)*(sigma/(b*omega)) * (2/(4-pi))^(1/3) / (gamma1^2)^(1/3) 
    tmp <- -sigma/(b*omega^2) * sign(gamma1)*(2*abs(gamma1)/(4-pi))^(1/3)
    dOmgdgam1 <- D23%*%Dtld33
    Adomgm1 <- array(dim=p)
    for (j in 1:p) 
      Adomgm1[j] <- tmp[j] * dOmgdgam1[ltdind(p)[j],j]/(2*omega[j])
    gamma1grad <- -2*(dAomgm1+Adomgm1)*OmgbIdelta + drop(vdOmgb%*%dOmgdgam1)    
    
    penaltygrad <- penaltygrad - PenF*(c2-c2tol)*c(mugrad,Sigmagrad,gamma1grad)
    if ( c2 < 0. || isTRUE(all.equal(c2,0.)) ) return(penaltygrad)     
  }  

  c1 <- sqrt((b2-(1.-b2)*beta02)/(1.+beta02))
  q1 <- 1./(c1*(1.+beta02))    
  q2 <- q1*(2*c1-q1)/2
  
  tmpM <- q1*q2*SigmaImu0%*%t(SigmaImu0)
  D32 <- -matrix( apply(kronecker(t(SigmaImu0),q1*SigmaI-tmpM),1,
                        rghtD,p=p,ncols=nvcovpar),
                  p,nvcovpar,byrow=TRUE ) 
  D33 <- q1*SigmaI-2*tmpM    
  eta <- q1*SigmaImu0
  list(ksi=ksi,Omega=Omega,eta=eta,OmegaInv=OmegaInv,
    nvcovpar=nvcovpar,penaltygrad=penaltygrad,
    D23=D23,D32=D32,D33=D33,Dtld32=Dtld32,Dtld33=Dtld33)
}
 
DPll.grad <- function(ksi,Omega,eta,OmegaInv,nvcovpar,penaltygrad,
  D23,D32,D33,Dtld32,Dtld33,
  y,n=ifelse(is.matrix(y),nrow(y),1),
  p=ifelse(is.matrix(y),ncol(y),length(y)),inoptm=FALSE,
#  PenF=0.,ldRtol=log(1E-5)+1-p,c2tol=1E-6)
  PenF=0.,ldRtol=log(1E-5)+1-p,c2tol=1E-6,limlnk2)
{
#  Thetagrad <- Theta.ll.grad(ksi,Omega,eta,y,n,p,OmegaInv)
  Thetagrad <- Theta.ll.grad(ksi,Omega,eta,y,n,p,OmegaInv,limlnk2=limlnk2)
  mugrad <- Thetagrad[1:p] 
  Omegagrad <- Thetagrad[(p+1):(p+nvcovpar)]
  etagrad <- Thetagrad[(p+nvcovpar+1):(2*p+nvcovpar)]
  
  Sigmagrad <- drop( -mugrad%*%Dtld32 + Omegagrad + Omegagrad%*%D23%*%Dtld32 +
                       etagrad%*%(D32+D33%*%Dtld32) )

  gamma1grad <- drop( -mugrad%*%Dtld33 + Omegagrad%*%D23%*%Dtld33 +
                        etagrad%*%D33%*%Dtld33 )
  penaltygrad + c(mugrad,Sigmagrad,gamma1grad)
}

OmgtoOmgbar.grad <- function(Omg,Omgdim,sclMat=NULL,Omgb=NULL) 
{
  # Note: Only three non-negative elemenst per row !!!
  #      Later try to find compact representations that multiply 
  #      these matrices more efficiently   
  
  omega <- sqrt(diag(Omg))
  if (is.null(sclMat)) sclMat <- 1./outer(omega,omega)
  if (is.null(Omgb)) Omgb <- sclMat * Omg
  
  Jacobdim <- Omgdim*(Omgdim+1)/2
  Jacob <- matrix(0.,Jacobdim,Jacobdim)
  diagind <- ltdind(Omgdim)
  Jind <- 0
  for (c in 1:(Omgdim-1)) {
    Jind <- Jind + 1 
    for (r in (c+1):Omgdim) {
      Jind <- Jind + 1 
      Jacob[Jind,Jind] <- sclMat[r,c]
      Jacob[Jind,diagind[r]] <- -Omgb[r,c]/(2*Omg[r,r])  
      Jacob[Jind,diagind[c]] <- -Omgb[r,c]/(2*Omg[c,c])
    }
  }
  Jacob
}  
