MINPHIARG <- -37	# minimum value allowed for the argument of a cumulative distribution of a standardized random variable 

## Things to do:  Check if these can now be replaced by the new functions of Azzalini's sn package!!

phioverPhi <- function(x,x0=MINPHIARG) if (x>x0) return(dnorm(x)/pnorm(x)) else return(dnorm(x0)/pnorm(x0))  
log2Phi <- function(x,x0=MINPHIARG) if (x>x0) return(log(2*pnorm(x))) else { Phi0 <- pnorm(x0) ; return(log(2*Phi0) + dnorm(x0)*(x-x0)/Phi0) }

ObsLik.ksigrad <- function(Xdevi,j,p,alpha,Omega)
{
   M <- matrix(0.,nrow=p,ncol=p)
   M[,j] <- -Xdevi
   phiarg <- (Xdevi/sqrt(diag(Omega)))%*%alpha
   -sum(diag(solve(Omega,(M+t(M)))) )/2 - phioverPhi(phiarg)*alpha[j]/sqrt(Omega[j,j]) 
}

ObsLik.alphagrad <- function(Xdevi,j,alpha,Omega)
{
   phiarg <- (Xdevi/sqrt(diag(Omega)))%*%alpha
   phioverPhi(phiarg)*Xdevi[j]/sqrt(Omega[j,j]) 
}

ObsLik.Omegagrad <- function(Xdevi,j1,j2,alpha,Omega)
{
   Vi <- outer(Xdevi,Xdevi)
   OmegaInv <- solve(Omega)
   grad <- -OmegaInv[j1,j2] + sum(diag(outer(OmegaInv[,j1],OmegaInv[,j2])%*%Vi))
   if (j1==j2)  {
   	phiarg <- (Xdevi/sqrt(diag(Omega)))%*%alpha
	grad <- grad/2 - phioverPhi(phiarg)*alpha[j1]*Xdevi[j1]/(2*Omega[j1,j1]^1.5) 
   }
   grad 
}
