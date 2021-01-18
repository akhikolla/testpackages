
fsf <- function(y, D, s = 20, K.max = 5, ddinv=NULL){
  m = dim(D)[1]
  if(is.null(ddinv)){
    ddt=D%*%t(D)
    if(rcond(ddt)<1e-12){
      ddinv <- Solve.banded(ddt+diag(rep(1e-6,m)),1,1, B = diag(1,dim(D)[1]))
    }else{
      ddinv <- Solve.banded(ddt,1,1, B = diag(1,dim(D)[1]))
    }
  }
  fit <- l0gen_c(y, D, s, K.max, ddinv)
  beta <- drop(fit$beta)
  z <- drop(fit$z)
  u <- drop(fit$u)
  return(list(y = y, beta = beta, v = z))
}

