fsfused.ada <- function(y, tau=1, s.max=20, T, eps=0.1){
  n <- length(y)
  m = n-1
  kk <- 1
  beta.all <- NULL
  while(kk*tau < s.max){

    s <- tau*kk
    re <- fsfused(y, s, T, K.max = 5)
    beta <- re$beta
    beta.all <- cbind(beta.all, beta)
    z <- re$z
    u <- re$u
    mse <- mean((y-beta)^2)
    #cat(T02, "\n")
    #cat(mse, "\n")
    if(mse < eps | s > s.max){
      break
    }else{
      kk <- kk + 1
    }
  }
  df <- 1:s
  return(list(y = y, beta = beta, v = z, beta.all=beta.all, df=df))
}





