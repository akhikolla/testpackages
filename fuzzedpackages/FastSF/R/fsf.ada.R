fsf.ada <- function(y, D, tau=1, s.max=20,  eps=0.1, ddinv=NULL){
  n <- length(y)
  m <- dim(D)[1]

  if(is.null(ddinv)){
    ddt=D%*%t(D)
    if(rcond(ddt)<1e-12){
      ddinv <- Solve.banded(ddt+diag(rep(1e-6,m)),1,1, B = diag(1,dim(D)[1]))
    }else{
      ddinv <- Solve.banded(ddt,1,1, B = diag(1,dim(D)[1]))
    }
  }

  kk <- 1
  beta.all <- NULL

  while(kk * tau < s.max){
    # cat(kk, "\n")
    s <- tau*kk
    re <- fsf(y = y, D = D, s = s, ddinv=ddinv)
    beta <- re$beta
    beta.all <- cbind(beta.all, beta)
    z <- re$z
    u <- re$u
    mse <- mean((y-beta)^2)
    # cat(mse, "\n")
    if(mse < eps | s > s.max){
      break
    }else{
      kk <- kk + 1
    }
  }
  df <- 1:s
  return(list(y = y, beta = beta, v = z, beta.all = beta.all, df = df))
}





