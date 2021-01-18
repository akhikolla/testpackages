ftf.ada <- function(y, k = 1, tau = 1, s.max=20, eps=0.1){
  n <- length(y)
  kk <- 1
  beta.all <- NULL
  while(kk * tau < s.max){
    # cat(kk, "\n")
    s <- kk*tau
    re <- ftf(y = y, k = k, s = s, K.max = 5)
    beta <- re$beta
    beta.all <- cbind(beta.all, beta)
    v <- re$v
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
  return(list(y = y, beta = beta, v = v, beta.all = beta.all, df = df))
}





