b2 <- 2./pi
ulbeta02 <- b2/(1.-b2) 

cfsbmu0 <- function(x,SigmaIchold,beta02tol=1E-3)
{
  p <- length(x)
  w <- array(dim=p)

  Resid <- ulbeta02 - beta02tol
  for (j in 1:p)  {
    wjcnt <- abs(x[j])*Resid
    w[j] <- sign(x[j])*sqrt(wjcnt)
    Resid <- Resid-wjcnt
  }
  backsolve(SigmaIchold,w)
}

getfskpar <- function(mu0,SigmaIchold,beta02tol=1E-3)
{
  p <- length(mu0)
  x <- array(dim=p)
  w <-  SigmaIchold %*% mu0
  Resid <- ulbeta02 - beta02tol
  for (j in 1:p)  {
    wj2 <- w[j]^2
    x[j] <- sign(w[j])*wj2/Resid
    Resid <- Resid-wj2
  }
  x  
}

