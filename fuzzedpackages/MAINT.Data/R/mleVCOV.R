#mleVCov <- function(mlegrad)  # Gets the Var-Cov of maximum likelihood estimators by the outer-gradient approximation
mleVCov <- function(mlegrad,limlnk2)  # Gets the Var-Cov of maximum likelihood estimators by the outer-gradient approximation

#  Arguments:

#  mlegrad   -  The gradients of each observation contribution to maximum likelihood function at its maximum
#		(a matrix with observations in rows and parameters in columns)

#  Value     - The asymptotic Var-Cov of the maximum likelihood estimators   
{
  if ( is.null(mlegrad) || any(!is.finite(mlegrad)) ) { return(NULL) } 
  npar <- ncol(mlegrad)	
  OtProd <- apply(mlegrad,1,function(v) outer(v,v))
  dim(OtProd) <- c(npar,npar,nrow(mlegrad))
  NegHessAp <- matrix(apply(OtProd,c(1,2),sum),nrow=npar,ncol=npar,byrow=FALSE)
#  pdwt.solve(NegHessAp,silent=TRUE)
  Safepdsolve(NegHessAp,maxlnk2=limlnk2,scale=TRUE)
}

NmleVCov <- function(nk,p,SigmaE,k=1,grpnames=NULL)  # Gets the Var-Cov of maximum likelihood estimators for Gaussian models

#  Arguments:

#  nk       -  Nunber of obervations (by group) used to get the mean and covariance matrix estimates
#  p        -  Nunber of variables
#  SigmaE   -  Matrix with covariance estimates
#  k        -  Number of different groups
#  grpnames -  Vector with the group names when k > 1

#  Value     - The Var-Cov matrix of the maximum likelihood estimators   
{
  if ( p!=nrow(SigmaE) || p!=ncol(SigmaE) ) {
    stop("NmleVCov arguments with wrong dimensions\n")
  }
  if (k!=length(nk)) {
    stop(paste("Dimension of argument n (",length(n),") does not agree ",
    "with the number of groups specified by argument k (",k,")\n",sep=""))
  }
  if ( k>1 && is.null(grpnames) ) {
    stop("grpnames argument is missing\n")
  }
  n <- sum(nk)
  npar <- k*p + p*(p+1)/2
  vnames <- dimnames(SigmaE)[[1]]
  mlevcov <- matrix(0.,nrow=npar,ncol=npar)	

      # Var-Covariance of sample means

  if (k==1)
  {
    parnames <- paste("mu_",vnames,sep="")
    mlevcov[1:p,1:p] <- SigmaE/n
  } else {
    parnames <- paste("mu_",vnames,"_",grpnames[1],sep="")
    mlevcov[1:p,1:p] <- SigmaE/nk[1]
    for (g in 2:k)
    {
      parnames <- c(parnames,paste("mu_",vnames,"_",grpnames[g],sep=""))
      muind <- (g-1)*p+1:p	
      mlevcov[muind,muind] <- SigmaE/nk[g]
    }
  }	

     # Var-Covariance of sample covariance matrix

  for (i in 1:p)
    parnames <- c(parnames,paste("Sigma_",vnames[i],"_",vnames[i:p],sep=""))
  ncnst <- (n-k)/n^2                    
  ind1 <- k*p
  for(i1 in 1:p) for(j1 in i1:p)  {
    ind1 <- ind1 + 1  
    ind2 <- k*p
    for(i2 in 1:p) for(j2 in i2:p)  {
      ind2 <- ind2 + 1
      if (ind2<=ind1)  {
        mlevcov[ind1,ind2] <- ncnst * ( SigmaE[i1,i2]*SigmaE[j1,j2] + SigmaE[i1,j2]*SigmaE[i2,j1] )
        if (ind2<ind1) mlevcov[ind2,ind1] <- mlevcov[ind1,ind2]
      }
    }
  }

  dimnames(mlevcov) <- list(parnames,parnames)
  mlevcov
}

HetSNmleVCov <- function(nk,p,k,Config,SigmaE,HomSNvcov,grpnames)  
# Gets the Var-Cov of maximum likelihood estimators for heterocedadtic SkewNormal models

#  Arguments:

#  nk          -  Nunber of obervations (by group) used to get the mean and covariance matrix estimates
#  p           -  Nunber of variables
#  k           -  Number of different groups
#  Config      -  Configuration of the Var-Covariance matrix
#  SigmaE      -  Matrix with covariance estimates
#  HomSNvcov   -  The parameter vcov matrix for the homoscedastic SkewNormal vector of residuals
#  grpnames    -  Vector with the group names

#  Value     - The Var-Cov matrix of the maximum likelihood estimators   
{
  if ( p!=nrow(SigmaE) || p!=ncol(SigmaE) ) {
    stop("HetSNmleVCov arguments with wrong dimensions\n")
  }
  if (k!=length(nk)) {
    stop(paste("Dimension of argument n (",length(n),") does not agree ",
    "with the number of groups specified by argument k (",k,")\n",sep=""))
  }
  n <- sum(nk)
  nmupar <- k*p
  SigConfind <- SigCind(Config,p/2)
  nSigpar <- length(SigConfind)
  Sigpar <- nmupar + 1:nSigpar
  gamma1par <- nmupar + nSigpar + 1:p 
  npar <- nmupar + nSigpar + p
  HmodSigpar <- p + 1:nSigpar 
  Hmodgamma1par <- p + nSigpar + 1:p 

  mlevcov <- matrix(0.,nrow=npar,ncol=npar)	
  vnames <- dimnames(SigmaE)[[1]]
  parnames <- paste("mu_",vnames,"_",grpnames[1],sep="")
  for (g in 2:k)  parnames <- c(parnames,paste("mu_",vnames,"_",grpnames[g],sep=""))
  for (i in 1:p)
    parnames <- c(parnames,paste("Sigma_",vnames[i],"_",vnames[i:p],sep=""))
  parnames <- parnames[c(1:nmupar,nmupar+SigConfind)]
  parnames <- c(parnames,paste("gamma1_",vnames,sep=""))
  dimnames(mlevcov) <- list(parnames,parnames)

  for (g1 in 1:k)
  {
    muind1 <- (g1-1)*p+1:p	
    mlevcov[muind1,muind1] <- SigmaE/nk[g1] + HomSNvcov[1:p,1:p]
    if (g1<k) for (g2 in (g1+1):k)
    {
      muind2 <- (g2-1)*p+1:p	
      mlevcov[muind1,muind2] <- mlevcov[muind2,muind1] <- HomSNvcov[1:p,1:p]
    }
    mlevcov[muind1,Sigpar] <- HomSNvcov[1:p,HmodSigpar]
    mlevcov[Sigpar,muind1] <- HomSNvcov[HmodSigpar,1:p]
    mlevcov[muind1,gamma1par] <- HomSNvcov[1:p,Hmodgamma1par]
    mlevcov[gamma1par,muind1] <- HomSNvcov[Hmodgamma1par,1:p]
  }	
  mlevcov[Sigpar,Sigpar] <- HomSNvcov[HmodSigpar,HmodSigpar]
  mlevcov[Sigpar,gamma1par] <- HomSNvcov[HmodSigpar,Hmodgamma1par]
  mlevcov[gamma1par,Sigpar] <- HomSNvcov[Hmodgamma1par,HmodSigpar]
  mlevcov[gamma1par,gamma1par] <- HomSNvcov[Hmodgamma1par,Hmodgamma1par]

  mlevcov
}

vcovCind <- function(Conf,p)
{
  if (!is.element(Conf,2:5)) {
    stop("Wrong value for Conf element\n")
  }
  q <- p/2
  mind <- NULL
  cind <- 0
  if (Conf==2)  {
    for (i in 1:q) {
      mind <- c(mind,cind+1:(q-i+1),cind+q+1)
      cind <- cind + q + (q-i+1)
    }
    for (i in 1:q) {
      mind <- c(mind,cind+1:(q-i+1))
      cind <- cind + (q-i+1)
    }
  }  
  if (Conf==3)  {
    for (i in 1:q) {
      mind <- c(mind,cind+1,cind+q+1)
      cind <- cind + q + (q-i+1)
    }
    for (i in 1:q) {
      mind <- c(mind,cind+1)
      cind <- cind + (q-i+1)
    }
  }  
  if (Conf==4)
  {
    for (i in 1:q) {
      mind <- c(mind,cind+1:(q-i+1))
      cind <- cind + q + (q-i+1)
    }
    for (i in 1:q) {
      mind <- c(mind,cind+1:(q-i+1))
      cind <- cind + (q-i+1)
    }
  }
  if (Conf==5) 
  {
    for (i in 1:q) {
      mind <- c(mind,cind+1)
      cind <- cind + q + (q-i+1)
    }
    for (i in 1:q) {
      mind <- c(mind,cind+1)
      cind <- cind + (q-i+1)
    }
  }
  mind  #  return(mind)
}

