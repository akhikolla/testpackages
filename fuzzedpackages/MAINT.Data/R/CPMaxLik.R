b2 <- 2./pi
ln2pi <- log(2*pi)

SNCMaxLik <- function(Data, Config, grouping=NULL, initpar=NULL, prevMod=NULL,
  EPS=1e-4,c2tol=1e-3,Beta02tol=1e-2,ldRtol=-500,limlnk2=log(1e8),OptCntrl=list())
{

#  Note: The Data argument should be a matrix containing the mid-points in the first columns and the log-ranges in the latter columns 

  if (!is.element(Config,2:5)) stop("Wrong value for Config argument\n")   
  maxsk <- 0.99527
  ulBeta02 <- b2/(1.-b2) - Beta02tol

  if (!is.matrix(Data)) Data <- as.matrix(Data)
  n <- nrow(Data)		    # number of observations
  p <- ncol(Data)		    # total number of variables (mid-points + log-ranges)
  q <- p/2			    # number of interval variables
  if (q==1) Config <- q1Config(Config) 
  if ( !is.null(grouping) )
  {                                                # Convert grouping factor into vector of integer indices with the second group
    grpind <- as.integer(grouping)-as.integer(2)   # indexed as 0 (for improved efficiency in C management of the beta2k matrix)
    k <- length(levels(grouping))                  # Number of groups
  } else {
    grpind <- rep(as.integer(-1),n)                # If the second group is indexed as 0, then the first group is indexed as -1
    k <- 1
  }
  nCoVarpar <- ncovp(Config,q,p)    # Number of free paramters in the VarCov matrix

  sd0_default <- 0.01
  if (!is.null(OptCntrl$sd0)) sd0 <- OptCntrl$sd0 
  else sd0 <- sd0_default 
  OptCntrl$sd0 <- NULL 
  parsd <- sd0*c(rep(1./sqrt(n),k*p),rep(0.1,nCoVarpar),rep(0.1,p))   # standard deviation hyper-parameters

  if (is.null(initpar))
  {
    if (is.null(prevMod)) prevMod <- SNCnf1MaxLik(Data,grouping,limlnk2=limlnk2,OptCntrl=OptCntrl)
    if (!is.null(prevMod$Sigma)) {
      Sigma <- CovC1toRestCov(p,q,prevMod$Sigma,Config)
    } else {
      Sigma <- diag(apply(Data,2,sd)^2)
    }  
    Spar <- GetCovPar(Sigma,Config=Config)
    if (is.null(Spar)) {
      cat("Data seems to be collinear (In GetCovPar 1: sample covariance matrix for Config",Config,"is not positive-definite)\n")
#      warning("Data seems to be collinear (sample covariance matrix is not positive-definite\n")
    } else {
      SigmaI <- Safepdsolve(Sigma,maxlnk2=limlnk2,scale=TRUE)
    } 
      if (is.null(Spar) || is.null(SigmaI)) {
      cat("Data seems to be collinear (In Safepdsolve: sample covariance matrix for Config",Config,"is not positive-definite)\n")
#      warning("Data seems to be collinear (sample covariance matrix is not positive-definite\n")
        return(list(lnLik=-Inf,ksi=NULL,Omega=NULL,Omega.cor=NULL,alpha=NULL,delta=NULL,
             mu=NULL,beta2k=NULL,Sigma=NULL,gamma1=NULL,admissible=NULL,c2=NULL,optres=NULL))
    }
    if (!is.null(prevMod$gamma1)) {
      gamma1 <- ifelse(abs(prevMod$gamma1)<maxsk,prevMod$gamma1,sign(prevMod$gamma1)*maxsk)
    } else {
      gamma1 <- rep(0.,p)
    }  
    mu0 <- sqrt(diag(Sigma)) * sign(gamma1)*(2*abs(gamma1)/(4.-pi))^(1/3) 
    Beta02 <- mu0 %*% drop(SigmaI %*% mu0)
    if (Beta02 > ulBeta02) gamma1 <- gamma1 * rep(sqrt(ulBeta02/Beta02)^3,length(gamma1)) 
    if ( is.null(grouping) ) {
      if (!is.null(prevMod$mu)) {
        prevModinitpar <- c(prevMod$mu,Spar,gamma1)
      } else {
        prevModinitpar <- c(colMeans(Data),Spar,gamma1)
      }  
    }  else  {
      if (!is.null(prevMod$mu)) {
        prevModinitpar <- c(prevMod$mu,prevMod$beta2k,Spar,gamma1)
      } else {
        prevModinitpar <- c(colMeans(Data),Spar,gamma1)   # Check later if this is indeed righ... is prevMod$mu really a vector?? 
      }  
    }

    prevModdev0 <- msnCP.dev(prevModinitpar,y=Data,grpind=grpind,Config=Config,k=k,limlnk2=limlnk2,nopenalty=TRUE)
    if (Config==2)  {
      C2res <- Cnf2MaxLik(Data,OptCntrl=OptCntrl)
      SigmaSrpar <- C2res$SigmaSr
    } else  {
      S <- var(Data)
      S[-RestCovInd(p,Config)] <- 0.
      SigmaSrpar <- GetCovPar(S,Config) 
      if (is.null(SigmaSrpar)) {
      cat("Data seems to be collinear (In GetCovPar 2: sample covariance matrix for Config",Config,"is not positive-definite)\n")
#        warning("Data seems to be collinear (Sample covariance matrix is not positive-definite)\n")
        return(list(lnLik=-Inf,ksi=NULL,Omega=NULL,Omega.cor=NULL,alpha=NULL,delta=NULL,
             mu=NULL,beta2k=NULL,Sigma=NULL,gamma1=NULL,admissible=NULL,c2=NULL,optres=NULL))
      }
    }
    if ( is.null(grouping) ) {
      NConfinitpar <- c(colMeans(Data),SigmaSrpar,rep(0.,p))
    }  else  {
      mug <- apply(Data,2,function(v) tapply(v,grouping,mean))
      mug1 <- mug[1,]
      if (k>2) {
        beta2k <- scale(mug[-1,],center=mug1,scale=FALSE)
      } else {
        beta2k <- mug[-1,] - mug1
      }
      NConfinitpar <- c(mug1,beta2k,SigmaSrpar,rep(0.,p))
    } 
    NConfdev0 <- msnCP.dev(NConfinitpar,y=Data,grpind=grpind,Config=Config,k=k,limlnk2=limlnk2,nopenalty=TRUE)
    if (NConfdev0 < prevModdev0) {
      initpar <- NConfinitpar 
      dev0 <- NConfdev0
    } else {
      initpar <- prevModinitpar 
      dev0 <- prevModdev0
    } 

  }  else  {
    dev0 <- msnCP.dev(initpar,y=Data,grpind=grpind,Config=Config,k=k,limlnk2=limlnk2,nopenalty=TRUE)
  }

  lb <- c(rep(-Inf,k*p),rep(EPS,p),rep(-Inf,nCoVarpar-p),rep(-maxsk,p))
  ub <- c(rep(Inf,k*p+nCoVarpar),rep(maxsk,p))

  if (is.null(OptCntrl$EnfCnstrs)) OptCntrl1 <- c(OptCntrl,EnfCnstrs=TRUE)
  else {
    OptCntrl1 <- OptCntrl
    OptCntrl1$EnfCnstrs <- TRUE
  }
  res <- RepLOptim(initpar,parsd,fr=msnCP.dev,gr=msnCP.dev.grad,lower=lb,upper=ub,
    y=Data,grpind=grpind,Config=Config,k=k,PenC=2*abs(dev0),limlnk2=limlnk2,c2tol=c2tol,ldRtol=ldRtol,control=OptCntrl1)

  if (dev0 < res$val) {
    res <- list(par=initpar,val=dev0,iterations=0,vallist=dev0,counts=0,convergence=NULL,message=NULL,
      hessian=NULL,hessegval=NULL,stderrors=NULL)
  }

  lnLik <- res$val/(-2)
  nSigmapar <- ncovp(Config,q,p)
  if (is.null(res$par)) {

   return( list(lnLik=lnLik,ksi=NULL,Omega=NULL,Omega.cor=NULL,alpha=NULL,delta=NULL,
      mu=NULL,Sigma=NULL,gamma1=NULL,admissible=FALSE,c2=NULL,optres=res)
    )
  }
  if ( is.null(grouping) ) {
    mu <- res$par[1:p]
    beta2k <- NULL
    Sigma <- RestCov(q,res$par[(p+1):(p+nSigmapar)],Config) 
    gamma1 <- res$par[(p+nSigmapar+1):(2*p+nSigmapar)]
    DP <- cnvCPtoDP(p,mu,Sigma,gamma1,limlnk2=limlnk2,silent=TRUE,tol=0.) 
    ksi <- DP$ksi
  }  else {
    mu <- matrix(nrow=k,ncol=p)
    mu[1,] <- mug1 <- res$par[1:p]
    beta2k <- matrix(res$par[(p+1):(k*p)],nrow=k-1,ncol=p)
#    mu[-1,] <- scale(beta2k,center=-mug1,scale=FALSE)
    if (k>2) {
      mu[-1,] <- scale(beta2k,center=-mug1,scale=FALSE)
    } else {
      mu[-1,] <- beta2k + mug1
    }
    Sigma <- RestCov(q,res$par[(k*p+1):(k*p+nSigmapar)],Config) 
    gamma1 <- res$par[(k*p+nSigmapar+1):((k+1)*p+nSigmapar)]
    DP <- cnvCPtoDP(p,mug1,Sigma,gamma1,limlnk2=limlnk2,silent=TRUE,tol=0.) 
    ksi <- matrix(nrow=k,ncol=p)
    ksi[1,] <- DP$ksi
    if (k>2) {
      ksi[-1,] <- scale(beta2k,center=-DP$ksi,scale=FALSE)
    } else {
      ksi[-1,] <- beta2k + DP$ksi
    }
  }  

  list(lnLik=lnLik,ksi=ksi,Omega=DP$Omega,Omega.cor=DP$Omega.cor,alpha=DP$alpha,delta=DP$delta,
    mu=mu,beta2k=beta2k,Sigma=Sigma,gamma1=gamma1,admissible=DP$admissible,c2=DP$c2,optres=res)
}

CovC1toRestCov <- function(p,q,Sigma,Config,rtol=1e-12)
{
  if (p!=2*q) stop("Argument p is not two times argument q\n")
  if (!is.element(Config,2:5)) stop("Wrong value for Config argument\n") 
  if (Config==5) return(diag(diag(Sigma)))
  if (Config==2)
  {
    if (q==1) return(Sigma)
    SigmaSr <- chol(Sigma)
    for (c in (q+1):p)
    {
      if (c!=1+q) SigmaSr[1,c] <- 0.
      for (r in 2:q) if (c!=r+q) 	     
        SigmaSr[r,c] <- -sum(SigmaSr[1:(r-1),r]*SigmaSr[1:(r-1),c])/SigmaSr[r,r]
    }
    RestSigma <- t(SigmaSr) %*% SigmaSr 
    atol <- rtol*max(abs(RestSigma))
    RestSigma[abs(RestSigma)<atol] <- 0.
    return(RestSigma)
  }
  RestSigma <- matrix(0.,p,p)
  if (Config==3) for (i in 1:q)
  {
    vind <- c(i,q+i)
    RestSigma[vind,vind] <- Sigma[vind,vind]
  }	
  if (Config==4)
  {
    RestSigma[1:q,1:q] <- Sigma[1:q,1:q]
    RestSigma[(q+1):p,(q+1):p] <- Sigma[(q+1):p,(q+1):p]
  }
  RestSigma  #  return(RestSigma)
}
