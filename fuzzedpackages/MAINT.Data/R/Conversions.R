b <- sqrt(2./pi)

cnvDPtoCP <- function(p,ksi,Omega,alpha)
{
   if (p==1)
   {
      omega <- sqrt(Omega) 
      delta <- alpha/sqrt(1.+alpha^2)
      muz <-  b*delta
      mu0 <- omega*muz
      mu <- ksi + mu0
      Sigma <- Omega - mu0^2
   }
   else  {
      omega <- sqrt(diag(Omega)) 
      Omegabar <- cov2cor(Omega)
      tmp <- drop(Omegabar %*% alpha)
      delta <- tmp/sqrt(rep(1.+drop(alpha%*%tmp),length(tmp)))
      muz <-  b*delta
      mu0 <- omega*muz

      mu <- ksi + mu0
      Sigma <- Omega - outer(mu0,mu0)
   }
   gamma1 <- ((4.-pi)/2) * (muz/sqrt(1.-muz^2))^3

   list(mu=mu,Sigma=Sigma,gamma1=gamma1)
}

cnvCPtoDP <- function(p, mu, Sigma, gamma1, limlnk2, silent=FALSE, tol=sqrt(.Machine$double.eps))
{
  c <- sign(gamma1)*(2*abs(gamma1)/(4.-pi))^(1./3)
  muz <- c/sqrt(1.+c^2)
  sigmaz <- sqrt(1. - muz^2)   
  delta <- muz/b
  if (p==1)
  {
    omega <- sqrt(Sigma)/sigmaz
    mu0 <- omega*muz

    ksi <- mu - mu0
    Omega <- Sigma + mu0^2
    Omegabar <- matrix(1.,1,1)
    c2 <- 1.-delta^2
    if (c2<0. || isTRUE(all.equal(c2,0.,tol=tol)))
      if (silent) return( 
        list(ksi=ksi,Omega=Omega,alpha=delta/tol,Omega.cor=Omegabar,delta=delta,c2=c2,
             admissible=FALSE,viol=tol-c2) )
    else stop("Inadmissible centred parameters\n")      
    alpha <- delta/sqrt(c2)
  }
  else  {
    omega <- sqrt(diag(Sigma))/sigmaz 
    mu0 <- omega*muz

    ksi <- mu - mu0
    Omega <- Sigma + outer(mu0,mu0)
    Omegabar <- cov2cor(Omega)
#    OmgbInv <- pdwt.solve(Omegabar,silent=TRUE)
#    OmgbInv <- Safepdsolve(Omegabar,maxlnk2=limlnk2,scale=FALSE)
    OmgbInv <- Safepdsolve(Omegabar,maxlnk2=limlnk2,scale=TRUE)
    if (is.null(OmgbInv)) {
      if (silent) return( 
        list(ksi=ksi,Omega=Omega,omega=omega,alpha=NULL,Omega.cor=Omegabar,delta=delta,c2=NULL,
             admissible=FALSE,viol=-determinant(Omegabar)$modulus) )
        else stop("Inadmissible centred parameters\n")       
    }   
    tmp <- drop(OmgbInv%*%delta)
    c2 <- 1. - drop(delta%*%tmp)
    if (c2<0. || isTRUE(all.equal(c2,0.,tol=tol)))
    {
      if (silent) {
        return( list(ksi=ksi,Omega=Omega,alpha=tmp/tol,Omega.cor=Omegabar,delta=delta,c2=c2,
          admissible=FALSE,viol=tol-c2) )
      } else {
        stop("Inadmissible centred parameters\n")
      }
    }      
    alpha <- tmp/rep(sqrt(c2),length(tmp))
  }
  
  list(ksi=ksi,Omega=Omega,omega=omega,alpha=alpha,Omega.cor=Omegabar,delta=delta,c2=c2,
       admissible=TRUE,viol=NULL)
}

