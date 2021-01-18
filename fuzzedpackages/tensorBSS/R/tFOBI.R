tFOBI <-
function(x, norm = NULL){
  
  if(is.null(norm)) norm <- rep(FALSE, length(dim(x)) - 1)
  
  if(length(dim(x)) == 2){
    returnlist <- FOBI(t(x))
    returnlist$S <- t(returnlist$S)
    returnlist2 <- list(S = returnlist$S, W = returnlist$W, Xmu = returnlist$Xmu, datatype = "iid")
    class(returnlist2) <- c("tbss", "bss") 
    return(returnlist2)
  }
  
  xmu <- apply(x, 1:(length(dim(x)) - 1), mean)
  
  stand <- tensorStandardize(x)
  x <- stand$x
  rotat <- tFOBIRotate(x, norm)
  x <- rotat$x
  
  r <- length(stand$S)
  W <- vector("list", r)
  for(i in 1:r){
    W[[i]] <- rotat$U[[i]]%*%stand$S[[i]]
  }
  
  returnlist <- list(S = x, W = W, norm = norm, Xmu = xmu, datatype = "iid")
  
  class(returnlist) <- c("tbss", "bss") 
  
  return(returnlist)
}
