tgFOBI <-
function(x, lags = 0:12, maxiter = 100, eps = 1e-06){
  
  xmu <- apply(x, 1:(length(dim(x)) - 1), mean)
  
  if(length(dim(x)) == 2){
    returnlist <- gFOBI(t(x), k = lags, eps = eps, maxiter = maxiter)
    returnlist$S <- t(returnlist$S)
    returnlist2 <- list(S = returnlist$S, W = returnlist$W, Xmu = xmu, datatype = "ts")
    class(returnlist2) <- c("tbss", "bss") 
    return(returnlist2)
  }
  
  stand <- tensorStandardize(x)
  x <- stand$x
  rotat <- tgFOBIRotate(x, lags, maxiter = maxiter, eps = eps)
  x <- rotat$x
  
  r <- length(stand$S)
  W <- vector("list", r)
  for(i in 1:r){
    W[[i]] <- rotat$U[[i]]%*%stand$S[[i]]
  }
  
  returnlist <- list(S = x, W = W, Xmu = xmu, datatype = "ts")
  
  class(returnlist) <- c("tbss", "bss") 
  
  return(returnlist)
}
