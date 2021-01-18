tJADE <-
function(x, maxiter = 100, eps = 1e-06){
  
  if(length(dim(x)) == 2){
    returnlist <- JADE(t(x), eps = eps, maxiter = maxiter)
    returnlist$S <- t(returnlist$S)
    returnlist2 <- list(S = returnlist$S, W = returnlist$W, Xmu = returnlist$Xmu, datatype = "iid")
    class(returnlist2) <- c("tbss", "bss") 
    return(returnlist2)
  }
  
  xmu <- apply(x, 1:(length(dim(x)) - 1), mean)
  
  stand <- tensorStandardize(x)
  x <- stand$x
  rotat <- tJADERotate(x, maxiter = maxiter, eps = eps)
  x <- rotat$x
  
  r <- length(stand$S)
  W <- vector("list", r)
  for(i in 1:r){
    W[[i]] <- rotat$U[[i]]%*%stand$S[[i]]
  }
  
  returnlist <- list(S = x, W = W, Xmu = xmu, datatype = "iid")
  
  class(returnlist) <- c("tbss", "bss") 
  
  return(returnlist)
}
