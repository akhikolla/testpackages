RestCov <- function(q,x,Config=1)
{
  if (Config==1) return(C1CovPar(q,x))
  if (Config==2) return(C2CovPar(q,x))
  if (Config==3) return(C3CovPar(q,x))
  if (Config==4) return(C4CovPar(q,x))
  if (Config==5) {
    if (length(x)!=2*q) 
      stop("C5CovPar: Wrong number of arguments for vector o Var-Cov parameters\n" )
    return(diag(x^2))
  }  
}

C1CovPar <- function(q,x) 
{
  p <- 2*q
  npar <- ncovp(1,q,p)
  if (length(x)!=npar) 
    stop("C1CovPar: Wrong number of arguments for vector o Var-Cov parameters\n" )
  A <- diag(x[1:p])
  A[col(A)>row(A)] <- x[(p+1):npar]
  t(A) %*% A
}   

C2CovPar <- function(q,x,rtol=1e-12) 
{
  p <- 2*q 
  npar <- ncovp(2,q,p)
  if (length(x)!=npar) 
    stop("C2CovPar: Wrong number of arguments for vector o Var-Cov parameters\n" )

  A <- diag(x[1:p])
  xind <- p
  for (c in 2:p) for(r in 1:(c-1))  {
    if ( (r<=q && c<=q)  || (r>q && c>q) || c==r+q )  
    {
      xind <- xind+1
      A[r,c] <- x[xind]
    }
    else if (r>1)   
      A[r,c] <- -sum(A[1:(r-1),r]*A[1:(r-1),c])/A[r,r]
  }  
  Sigma <- t(A) %*% A
  atol <- rtol*max(abs(Sigma))
  Sigma[abs(Sigma)<atol] <- 0.
  Sigma
} 
   
C3CovPar <- function(q,x,rtol=1e-12) 
{
  p <- 2*q 
  npar <- ncovp(3,q,p)
  if (length(x)!=npar) 
    stop("C3CovPar: Wrong number of arguments for vector o Var-Cov parameters\n" )
  Sigma <- matrix(0.,nrow=p,ncol=p)
  A <- matrix(0.,nrow=2,ncol=2)
  for (i in 1:q) {
     A[1,1] <- x[i]
     A[2,2] <- x[q+i]
     A[1,2] <- x[p+i]
     iind <- c(i,q+i)
     Sigma[iind,iind] <- t(A)%*%A
  }
  atol <- rtol*max(abs(Sigma))
  Sigma[abs(Sigma)<atol] <- 0.
  Sigma
}

C4CovPar <- function(q,x,rtol=1e-12) 
{
  p <- 2*q 
  npar <- ncovp(4,q,p)
  if (length(x)!=npar) 
    stop("C4CovPar: Wrong number of arguments for vector o Var-Cov parameters\n" )
  if (q==1) return(diag(x^2))
  
  RCCind <- (p+1):(p+q*(q-1)/2)  
  RLrLrind <- (p+q*(q-1)/2+1):(p+q*(q-1))  
  Sigma <- matrix(0.,nrow=p,ncol=p)
  A <- diag(x[1:q])
  A[col(A)>row(A)] <- x[RCCind]
  Sigma[1:q,1:q] <- t(A)%*%A
  diag(A) <- x[(q+1):(2*q)]
  A[col(A)>row(A)] <- x[RLrLrind]
  Sigma[(q+1):p,(q+1):p] <- t(A)%*%A
  atol <- rtol*max(abs(Sigma))
  Sigma[abs(Sigma)<atol] <- 0.
  Sigma
}   

RestCovInd <- function(p,Config)
{
  if (Config==1) return(1:p^2)
  if (Config==2) return(C2RestCovInd(p))
  if (Config==3) return(C3RestCovInd(p))
  if (Config==4) return(C4RestCovInd(p))
  if (Config==5) return(C5RestCovInd(p))
}

C2RestCovInd <- function(p) 
{
  q <- p/2 
  ind <- NULL
  for (c in 1:q) ind <- c(ind,(c-1)*p+c((1:q),q+c))
  for (c in (q+1):p) ind <- c(ind,(c-1)*p+c(c-q,q+(1:q)))
  ind
} 

C3RestCovInd <- function(p) 
{
  q <- p/2   
  ind <- NULL
  for (c in 1:q) ind <- c(ind,(c-1)*p+c(0,q)+c)
  for (c in (q+1):p) ind <- c(ind,(c-1)*p+c(-q,0)+c)
  ind
}   
   
C4RestCovInd <- function(p) 
{
  q <- p/2 
  ind <- NULL
  for (c in 1:q) ind <- c(ind,(c-1)*p+1:q)
  for (c in (q+1):p) ind <- c(ind,(c-1)*p+(q+1):p)
  ind
}

C5RestCovInd <- function(p) (0:(p-1))*p + 1:p

ncovp <- function(Config,q,p=2*p)
{
  if (Config==1) return(p*(p+1)/2)
  if (Config==2) return(p+q+q*(q-1))
  if (Config==3) return(p+q)
  if (Config==4) return(p+q*(q-1))
  if (Config==5) return(p)                                              
}


GetCovPar <- function(S,Config=1,test=TRUE) 
{
  if (!is.matrix(S) || !is.element(Config,1:5))  stop("Wrong argument types\n")
  p <- nrow(S)
  q <- p/2
  if (ncol(S)!=p) stop("Matrix S is not squared\n")
  
  if (Config==1 || Config==2) S1 <- S
  if (Config==3 || Config==5) S1 <- diag(diag(S))
  if (Config==3) for (v in 1:q) S1[v,q+v] <- S1[q+v,v] <- S[q+v,v]
#  if (Config==4) S1 <- rbind(cbind(S[1:q,1:q],matrix(0.,q,q)),cbind(matrix(0.,q,q),S[(q+1):p,(q+1):p]))
  if (Config==4) {
    S1 <- rbind.data.frame(
      cbind.data.frame(S[1:q,1:q],matrix(0.,q,q,dimnames=list(rownames(S)[1:q],colnames(S)[(q+1):p]))),
      cbind.data.frame(matrix(0.,q,q,dimnames=list(rownames(S)[(q+1):p],colnames(S)[1:q])),S[(q+1):p,(q+1):p])
    )
  }
    
  if (test) {
    Sr <- try(chol(S1, pivot = FALSE), silent = TRUE)
    if(class(Sr)[1] == "try-error") return(NULL)
  } else {
    Sr <- chol(S)
  }
  if (Config==1) return(c(diag(Sr),Sr[col(Sr)>row(Sr)]))
  if (Config==5) return(diag(Sr))
  
  q <- p/2
  npar <- ncovp(Config,q,p)
  retpar <- array(dim=npar)  
  retpar[1:p] <- diag(Sr)  
  rtpind <- p
  for (c in 2:p) for(r in 1:(c-1))  {
    if (  ( c==r+q && (Config==2 || Config==3) ) ||  
          ( ((r<=q && c<=q) || (r>q && c>q)) && (Config==2 || Config==4) ) 
    )  
    {
      rtpind <- rtpind+1
      retpar[rtpind] <- Sr[r,c]      
    }
  }
  retpar
} 

GetSigmaPar <- function(S,Config=1) 
{
  if (!is.matrix(S) || !is.element(Config,1:5))  stop("Wrong argument types\n")
  p <- nrow(S)
  q <- p/2
  if (ncol(S)!=p) stop("Matrix S is not squared\n")
  
  if (Config==1 || Config==2) S1 <- S
  if (Config==3 || Config==5) S1 <- diag(diag(S))
  if (Config==3) for (v in 1:q) S1[v,q+v] <- S1[q+v,v] <- S[q+v,v]
  if (Config==4) {
    S1 <- rbind.data.frame(
      cbind.data.frame(S[1:q,1:q],matrix(0.,q,q,dimnames=list(rownames(S)[1:q],colnames(S)[(q+1):p]))),
      cbind.data.frame(matrix(0.,q,q,dimnames=list(rownames(S)[(q+1):p],colnames(S)[1:q])),S[(q+1):p,(q+1):p])
    )
  }
    
  if (Config==1) return(c(diag(S1),S1[col(S1)>row(S1)]))
  if (Config==5) return(diag(S1))
  
  q <- p/2
  npar <- ncovp(Config,q,p)
  retpar <- array(dim=npar)  
  retpar[1:p] <- diag(S1)  
  rtpind <- p
  for (c in 2:p) for(r in 1:(c-1))  {
    if (  ( c==r+q && (Config==2 || Config==3) ) ||  
          ( ((r<=q && c<=q) || (r>q && c>q)) && (Config==2 || Config==4) ) 
    )  
    {
      rtpind <- rtpind+1
      retpar[rtpind] <- S1[r,c]      
    }
  }
  retpar
} 


ltind1 <- function(p,r,c) { # ordering lower-trg vec indices including diagonal 
  if (r < c || r > p) stop("Wrong argument values for ltind1\n")
  d <- p-c+1 
  (p*(p+1)-d*(d+1))/2 + r-c+1 
}

ltind0 <- function(p,r,c) { # ordering lower-trg vec indices excluding diagonal 
  if (r <= c || r > p) stop("Wrong argument values for ltind0\n")
  d <- p-c+1 
  (p*(p-1)-d*(d-1))/2 + r-c 
}

utind1 <- function(p,r,c) { # ordering upper-trg vec indices including diagonal 
  if (r > c || c > p) stop("Wrong argument values for utind1\n")
  (c-1)*c/2 + r 
}

utind0 <- function(p,r,c) { # ordering upper-trg vec indices excluding diagonal 
  if (r >= c || c > p) stop("Wrong argument values for utind0\n")
  (c-2)*(c-1)/2 + r 
}

lttoutind <- function(p) 
{  
  # converts indices of stacked vectors of different elements of 
  # symmetric matrices from lower-triangular to upper-triangular 
  
  res <- array(dim=p*(p+1)/2)
  ind <- 0
  for (c in 1:p) for (r in c:p) {
    ind <- ind +1
    res[ind] <- utind1(p,c,r)
  }
  res
}  

uttoltind <- function(p) 
{  
  # converts indices of stacked vectors of different elements of 
  # symmetric matrices from upper-triangular to lower-triangular 
  
  res <- array(dim=p*(p+1)/2)
  ind <- 0
  for (c in 1:p) for (r in 1:c) {
    ind <- ind +1
    res[ind] <- ltind1(p,c,r)
  }
  res
}  


c2cind <- function(q,p) 
{  
  if (p!=2*q) stop("Wrong arguments for c2cind function\n")
  if (q==1) return(1:2)
  qtrdim <- q*(q-1)/2
  ind <- c(1:p,p+(1:qtrdim))
  cind <- p+qtrdim
  for (c in (q+1):p) for (r in 1:(c-1))  {
    cind <- cind + 1
    if (c==q+r || r>q) ind <- c(ind,cind)
  }
  ind
}  

c2rfind <- function(q,p) 
{  
  if (p!=2*q) stop("Wrong arguments for c2rfind function\n")
  if (q==1) return(0)
  ind <- NULL
  cind <- q*(q+1)/2
  for (c in (q+1):p) {
    for (r in 1:q)  if (c!=q+r) ind <- c(ind,cind+r)
    cind <- cind + c
  }
  ind
}  

RestCov.grad <- function(q,x,Config)
{  
  p <- 2*q
  npart <- p*(p+1)/2
  pparp <- ncovp(Config,q,p)
  
  if (Config==1) Jacob <- C1CPgrad(x,p,npart,pparp)
  if (Config==2) Jacob <- C2CPgrad(x,q,p,npart)
  if (Config==3) Jacob <- C3CPgrad(x,q,p,npart,pparp)
  if (Config==4) Jacob <- C4CPgrad(x,q,p,npart,pparp)
  if (Config==5) Jacob <- C5CPgrad(x,p,npart,pparp)
  
  Jacob
}

C1CPgrad <- function(x,p,npart,pparp)
{
  A <- diag(x[1:p])
  A[row(A)<col(A)] <- x[(p+1):pparp]
  Jacob <- matrix(0.,npart,pparp)
  
  ind <- 0
  for (c in 1:p)  {
    if (c>1) for (r in 1:(c-1))  {
      ind <- ind+1
      for (l in 1:r)  {
        Jacob[ind,p+utind0(p,l,c)] <- A[l,r]
        if (l<r) Jacob[ind,p+utind0(p,l,r)] <- A[l,c]
        else Jacob[ind,r] <- A[l,c]
      }  
    }
    ind <- ind+1
    Jacob[ind,c] <- 2*A[c,c]
    if (c>1) for (l in 1:(c-1))  Jacob[ind,p+utind0(p,l,c)] <- 2*A[l,c]
  }  
  
  Jacob  
}

C2CPgrad <- function(x,q,p,npart) 
{
  Sigma <- RestCov(q,x,Config=2)
  SgmSr <- chol(Sigma)  
  tmp <- SgmSr[row(SgmSr)<=col(SgmSr)]
  dind <- utdind(p)
  fullpar <- c(tmp[dind],tmp[-dind])
  A <- diag(fullpar[1:p])
  A[row(A)<col(A)] <- fullpar[(p+1):npart]
  Jacob <- C1CPgrad(fullpar,p,npart,npart)
  if (q==1) return(Jacob)
  
  for (r in (q-1):1)  {
    for (c in (q+1):p) for (i in (r+1):q) if (c!=q+i)  {
      cind <- p+utind0(p,r,c)
      Jacob[,cind] <- 
        Jacob[,cind] - (A[r,i]*Jacob[,p+utind0(p,i,c)])/A[i,i]
    }  
    for (c in (r+1):q) for (j in (q+1):p) if (j!=q+c)  {
      cind <- p+utind0(p,r,c)
      Jacob[,cind] <- 
        Jacob[,cind] - (A[r,j]*Jacob[,p+utind0(p,c,j)])/A[c,c]
    }
  }  
  for (c in q:2)  
    for (j in (q+1):p) if (j!=q+c)
      Jacob[,c] <-
    Jacob[,c] + sum(A[1:(c-1),c]*A[1:(c-1),j])*Jacob[,p+utind0(p,c,j)]/A[c,c]^2 
  
  Jacob[c2rfind(q,p),] <- 0.
  Jacob[,c2cind(q,p)]
}  

C3CPgrad <- function(x,q,p,npart,pparp)
{
  Jacob <- matrix(0.,npart,pparp)
  dind <- utdind(p)
  for (v in 1:q) {
    Jacob[dind[v],v] <- 2*x[v]
    Jacob[dind[q+v],q+v] <- 2*x[q+v]
    v12rind <- utind1(p,v,v+q) 
    v12cind <- p + v 
    Jacob[v12rind,v] <- x[v12cind]
    Jacob[v12rind,v12cind] <- x[v]
    Jacob[dind[q+v],v12cind] <- 2*x[v12cind]
  }  
  Jacob  
}

C4CPgrad <- function(x,q,p,npart,pparp)
{
  if (q==1) return(matrix(c(2*x[1],rep(0.,4),2*x[2]),3,2))
  
  MidPind <- c(1:q,(p+1):(p+q*(q-1)/2))
  LRind <- c((q+1):p,(p+q*(q-1)/2+1):(q*(q+1)))
  qnpart <- q*(q+1)/2
  JacobMidP <- C1CPgrad(x[MidPind],q,qnpart,qnpart)
  JacobLR <- C1CPgrad(x[LRind],q,qnpart,qnpart)
  
  Jacob <- matrix(0.,npart,pparp)  
  Jacob[1:qnpart,MidPind] <- JacobMidP
  frind <- qnpart
  prind <- 0
  for (c in 1:q) {
    frind <- frind + c-1 + q
    prind <- prind + c-1 
    for (r in 1:c)  Jacob[frind+r,LRind] <- JacobLR[prind+r,]
  }  
  Jacob  
}

C5CPgrad <- function(x,p,npart,pparp)
{
  Jacob <- matrix(0.,npart,pparp)
  Jind <- utdind(p)
  for (i in 1:p) Jacob[Jind[i],i] <- 2*x[i]
  
  Jacob  
}

RestCovj <- function(x,q,Config=1,j)  
{  
  Mtmp <- RestCov(q,x,Config)
  as.vector(Mtmp[col(Mtmp)>=row(Mtmp)])[j]
}


