#MStep <- function(X,z,Cf,Homoc,k2max,tautol,invalcode=NULL)
MStep <- function(X,z,Cf,Homoc,k2max,tautol,MaxVarGRt,invalcode=NULL)
{
   n <- nrow(X)
   p <- ncol(X)
   q <- p/2
   k <- ncol(z)

   midpind <- 1:q 
   lrind <- (q+1):p
   Obsnames <- rownames(X)
   Vnames <- colnames(X)
   grpnames <- colnames(z)

   nktol <- tautol/n
   nk <- apply(z,2,sum)
   if ( any(!is.finite(nk)) || any(nk < nktol) )  return(invalcode) 	
   tau <- nk/n
   muk <- matrix(nrow=p,ncol=k,dimnames=list(Vnames,grpnames))	
   Wk <- array(0.,dim=c(p,p,k),dimnames=list(Vnames,Vnames,grpnames))
   if (Homoc==FALSE) Sigmak <- array(dim=c(p,p,k),dimnames=list(Vnames,Vnames,grpnames))
   zsqrt <- sqrt(z)

   for (g in 1:k)  {
     muk[,g] <- apply(sweep(X,1,z[,g],FUN="*"),2,sum)/nk[g]
     wdev <-  sweep(scale(X,center=muk[,g],scale=FALSE),1,zsqrt[,g],FUN="*")
     if (Cf==1) Wk[,,g] <- t(wdev)%*%wdev
     if (Cf==3)  {
       Wk[midpind,midpind,g] <- t(wdev[,midpind]) %*% wdev[,midpind] 
       Wk[lrind,lrind,g] <- t(wdev[,lrind]) %*% wdev[,lrind] 
     }
     if (Cf==2 || Cf==4 ) Wk[,,g] <- diag(apply(wdev^2,2,sum))
     if (Cf==2) for (j in 1:q) Wk[j,q+j,g] <- Wk[q+j,j,g] <- sum(wdev[,j]*wdev[,q+j])
     if (Homoc==FALSE) Sigmak[,,g] <- Wk[,,g]/nk[g]
  }

  if (Homoc==TRUE) {
    Sigma <- apply(Wk,c(1,2),sum)/n
    if (CheckSigmaSing(Cf,Sigma,limlnk2=log(k2max),scale=TRUE)) return(invalcode)  
     dimnames(Sigma) <- list(Vnames,Vnames)
    return(list(tau=tau,muk=muk,Sigma=Sigma,Sigmak=NULL))
  } else {  
#    if (CheckSigmak(Cf,Sigmak,MaxVarGRt=MaxVarGRt,limlnk2=log(k2max),scale=TRUE)) return(invalcode)
    if ( CheckSigmak(Cf,Sigmak,MaxVarGRt=MaxVarGRt,limlnk2=log(k2max),scale=TRUE)[[1]] != "Valid" ) return(invalcode)
    return(list(tau=tau,muk=muk,Sigma=NULL,Sigmak=Sigmak))
  }
}

