Pertubtau <- function(tau,pertubtau,EPS=1e-6)
{
  halfp <- pertubtau/2
  tau <- pmax(tau + runif(length(tau),-halfp,halfp),EPS)
  tau/sum(tau)
}  

Pertubz <- function(z,pertubz,EPS=1e-6)
{
   for (g in 1:ncol(z)) {
     halfp <- pertubz[,g]/2
     z[,g] <- pmax(z[,g] + runif(nrow(z),-halfp,halfp),EPS)
   }
   zscaled <- t(scale(t(z),center=FALSE,scale=rowSums(z)))
   attr(zscaled,"scaled:scale") <- NULL
   zscaled
}  

RepEMGauss <- function(X, n, p, k, ISol, pertub, nrep, Cf, Homoc,
  maxiter, convtol, protol, k2max, MaxSctEgvlRt, SctEgvCnstr, MaxVarGRt,
  Ifact=1e-3, Vfact=1e-1, maxinittrials=20)
{
  INVALID <- -1e99
  Lik <- numeric(nrep)

  Pertubate <- function()
  {
    Rdiagpert <- 1e-6 
    epsilon <- Rdiagpert^2

    Vnames <- names(X)
    Onames <- rownames(X)
    inittrial <- 0
    LnLik <- INVALID
    q <- p/2
    if (Homoc==TRUE) Sigmak <-NULL
    else {  
      Sigma <- NULL
      Sigmak <- array(dim=c(p,p,k),dimnames=list(Vnames,Vnames,NULL))
    }  
    if (!is.null(pertub$tau)) tau <- numeric(k)
    z <- ISol$z

    while (LnLik==INVALID && inittrial<maxinittrials)
    {
      if (!is.null(pertub$z))  {
        z <- Pertubz(z,pertub$z)
        LnLik <- MClusLikz(z,X,n,p,k,Homoc,k2max=k2max,penalty=INVALID) 
      } else {
        if (!is.null(pertub$tau)) {
          tau <- Pertubtau(ISol$tau,pertub$tau)
        } else {
          tau <- ISol$tau
        }  
        if (!is.null(pertub$muk)) {
          muk <- matrix(rnorm(k*p,ISol$muk,pertub$muk),nrow=p,ncol=k,dimnames=list(Vnames,NULL))
        } else {
          muk <- ISol$muk
        } 
      
        if (!is.null(pertub$Stdev)) {
          q0 <- p*(p-1)/2
       	  if (Homoc==TRUE) {
            if (CheckSigmaSing(Cf,ISol$Sigma,limlnk2=log(k2max),scale=TRUE)==TRUE) {
              diag(ISol$Sigma) <- (1.+Ifact)*diag(ISol$Sigma)
            }
       	    Sk <- 1
       	  } else {
            SigmaType <- list("Unknown")
            while (SigmaType[[1]] != "Valid") {
              SigmaType <- CheckSigmak(Cf,ISol$Sigmak,limlnk2=log(k2max),MaxVarGRt=MaxVarGRt,scale=TRUE)             
              if (SigmaType[[1]]=="Singular") {
                for (g in 1:k) diag(ISol$Sigmak[,,k]) <- (1+Ifact)*diag(ISol$Sigmak[,,k])
              } else if (SigmaType[[1]]=="LargeVarR") {  # Too large across groups variance ratio for the variables in SigmaType$viollst
                cnstf <- (1+Vfact)*SigmaType$Maxviol
                for (v in SigmaType$viollst) {
                  mingv <- SigmaType$mingrpVar[v]
                  ISol$Sigmak[v,v,mingv] <- cnstf*ISol$Sigmak[v,v,mingv]
                } 
              }
            }  
       	    Sk <- k
       	  }  
          if (Cf==2)  RSr <- vector("list",q)
          else if (Cf==3) RSr <- vector("list",2) 
          for (g in 1:Sk)  {
            if (Homoc==TRUE) {
              Stdev <- sqrt(diag(ISol$Sigma))
              invalidStdev <- which(!is.finite(Stdev) | Stdev <= 0.)  
              Stdev[invalidStdev] <- epsilon
              pertSigma <- ISol$Sigma + Rdiagpert*Stdev*diag(1.,p)
              if (Cf==1) RSr <- safecholR(pertSigma)
              else if (Cf==2) for (v in 1:q) RSr[[v]] <- safecholR(pertSigma[c(v,q+v),c(v,q+v)])
              else if (Cf==3) { 
                RSr[[1]] <- safecholR(pertSigma[1:q,1:q])
                RSr[[2]] <- safecholR(pertSigma[q+1:q,q+1:q])
              }
            } else {
              Stdev <- sqrt(diag(ISol$Sigmak[,,g]))
              invalidStdev <- which(!is.finite(Stdev) | Stdev <= 0.)  
              Stdev[invalidStdev] <- epsilon
              pertSigmakg <- ISol$Sigmak[,,g] + Rdiagpert*Stdev*diag(1.,p)
              if (Cf==1) RSr <- safecholR(pertSigmakg)
              else if (Cf==2) for (v in 1:q) RSr[[v]] <- safecholR(pertSigmakg[c(v,q+v),c(v,q+v)])
              else if (Cf==3) { 
                RSr[[1]] <- safecholR(pertSigmakg[1:q,1:q])
                RSr[[2]] <- safecholR(pertSigmakg[q+1:q,q+1:q])
              }
            }  
            Stdev <- abs(rnorm(p,Stdev,pertub$Stdev[(g-1)*p+1:p]))
            if (Cf==4) {
              if (Homoc==TRUE) Sigma <- diag(Stdev^2)
              else Sigmak[,,g] <- diag(Stdev^2)
            } else {
              if (Cf==1) {
                uptrigel <- upper.tri(RSr) 
                RSr[uptrigel] <- rnorm(q0,RSr[uptrigel],pertub$cor[(g-1)*q0+1:q0])
                R <- cov2cor(t(RSr) %*% RSr)        
                if (Homoc==TRUE) {
                  Sigma <- R * outer(Stdev,Stdev)
                } else {
                  Sigmak[,,g] <- R * outer(Stdev,Stdev)
                }
              } else if (Cf==2 || Cf==3) {            
                if (Homoc==TRUE) Sigma <- matrix(0.,nrow=p,ncol=p)
                else Sigmak[,,g] <- 0.
                if (Cf==2) {            
                  for (v in 1:q) {            
                    RSr[[v]][1,2] <- rnorm(1,RSr[[v]][1,2],pertub$cor[(g-1)*q+v])
                    R <- cov2cor(t(RSr[[v]]) %*% RSr[[v]])
                    indices <- c(v,q+v)         
                    if (Homoc==TRUE) {
                      Sigma[indices,indices] <- R * outer(Stdev[indices],Stdev[indices])
                    } else {
                      Sigmak[indices,indices,g] <- R * outer(Stdev[indices],Stdev[indices])
                    }
                  }
                } else if (Cf==3) {     
                  q4 <- q*(q-1)/2       
                  uptrigel <- upper.tri(RSr[[1]])  # it could just as well be upper.tri(RSr[[2]])... 
                  for (b in 1:2) {            
                    RSr[[b]][uptrigel] <- rnorm(q4,RSr[[b]][uptrigel],pertub$cor[(g-1)*2*q4+(b-1)*q4+1:q4])
                    R <- cov2cor(t(RSr[[b]]) %*% RSr[[b]])
                    indices <- (b-1)*q+1:q         
                    if (Homoc==TRUE) {
                      Sigma[indices,indices] <- R * outer(Stdev[indices],Stdev[indices])
                    } else {
                      Sigmak[indices,indices,g] <- R * outer(Stdev[indices],Stdev[indices])
                    }
                  }
                }
              }  
            }
          }
          if (Homoc==TRUE) dimnames(Sigma) <- list(Vnames,Vnames) 

        }  else { 
          Sigma <- ISol$Sigma ; Sigmak <- ISol$Sigmak
        } 
        if (Homoc==TRUE) {
          LnLik <- MClusLikpar(X,n,p,k,tau=tau,muk=muk,Homoc=TRUE,Sigma=Sigma,k2max=k2max,penalty=INVALID)
        } else {
          LnLik <- MClusLikpar(X,n,p,k,tau,Homoc=FALSE,muk=muk,Sigmak=Sigmak,k2max=k2max,penalty=INVALID)
        }
      } 
      inittrial <- inittrial+1  
    }
    return( list(z=z,tau=tau,muk=muk,Sigma=Sigma,Sigmak=Sigmak,LnLik=LnLik) ) 
  }

  Replicate <- function(z,tau,muk,Sigma,Sigmak,LnLik)
  {
    Vnames <- names(X)
    Onames <- rownames(X)
    Cnames <- paste("CP",1:k,sep="")

    if (Homoc==TRUE) {
              
      HomSol <- list(tau=tau,muk=muk,Sigma=Sigma,Sigmak=NULL,z=z,clusters=NULL,LnLik=LnLik,npar=NULL,BIC=NULL,AIC=NULL )
      res <- .Call( "CEMGauss", as.matrix(X), k, Cf, TRUE, 
        maxiter, protol, convtol, k2max, MaxSctEgvlRt, 
        z,tau,muk,Sigma,NULL,LnLik, FALSE, SctEgvCnstr, MaxVarGRt, 
        PACKAGE = "MAINT.Data" )
      if ( is.null(res$clusters) || any(!is.finite(res$tau)) || any(!is.finite(res$muk)) || any(!is.finite(res$Sigma)) || 
           CheckSigmaSing(Cf,res$Sigma,limlnk2=log(k2max),scale=TRUE)==TRUE )
      {        
        res$LnLik <- -Inf        
      }  else {
        rownames(res$muk) <- rownames(res$Sigma) <- colnames(res$Sigma) <- Vnames
        names(res$clusters) <- Onames
        rownames(res$z) <- Onames
        names(res$tau) <- colnames(res$muk) <- colnames(res$z) <- Cnames
        res$clusters <- Cnames[res$clusters]
      }
      
    } else {
       
      HetSol <- list(tau=tau,muk=muk,Sigma=NULL,Sigmak=Sigmak,z=z,clusters=NULL,LnLik=LnLik,npar=NULL,BIC=NULL,AIC=NULL )
      res <- .Call( "CEMGauss", as.matrix(X), k, Cf, FALSE, 
        maxiter, protol, convtol, k2max, MaxSctEgvlRt, 
        z,tau,muk,NULL,Sigmak,LnLik, FALSE, SctEgvCnstr, MaxVarGRt,
        PACKAGE = "MAINT.Data" )
      if ( is.null(res$clusters) || any(!is.finite(res$tau)) || any(!is.finite(res$muk)) || any(!is.finite(res$Sigmak)) || 
           CheckSigmak(Cf,res$Sigmak,limlnk2=log(k2max),MaxVarGRt=MaxVarGRt,scale=TRUE)[[1]] != "Valid" )
      {
        res$LnLik <- -Inf  
      }  else {
        rownames(res$muk) <- dimnames(res$Sigmak)[[1]]  <- dimnames(res$Sigmak)[[2]] <- Vnames
        names(res$clusters) <- Onames
        dimnames(res$Sigmak)[[3]] <- Cnames
        rownames(res$z) <- Onames
        names(res$tau) <- colnames(res$muk) <- colnames(res$z) <- Cnames
        res$clusters <- Cnames[res$clusters]
      }       
    }
    res # return(res)
  }

  BestLnLik <- ISol$LnLik
  for (rep in 1:nrep)  {
    PertSol <- Pertubate()
    if (PertSol$LnLik==INVALID) {
      Lik[rep] <- -Inf
    } else { 
      NxtSol <- Replicate(PertSol$z,PertSol$tau,PertSol$muk,PertSol$Sigma,PertSol$Sigmak,PertSol$LnLik)
      Lik[rep] <- NxtSol$LnLik
      if (!is.null(NxtSol$LnLik) && is.finite(NxtSol$LnLik) && NxtSol$LnLik > BestLnLik) {
        ISol <- NxtSol
        BestLnLik <- NxtSol$LnLik
      }
    }     
  }

  return(list(BestSol=ISol,alllnLik=Lik))
} 

