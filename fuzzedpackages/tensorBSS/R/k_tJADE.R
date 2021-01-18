k_tJADE <-
  function(x, k = NULL, maxiter = 100, eps = 1e-06){
    
    if(is.null(k)){
      k <- rep(1, length(dim(x)) - 1)
    }
    
    r <- length(dim(x)) - 1
    
    if(r == 1){
      returnlist <- k_JADE(t(x), k = k[1], eps = eps, maxiter = maxiter)
      returnlist$S <- t(returnlist$S)
      returnlist2 <- list(S = returnlist$S, W = returnlist$W, Xmu = returnlist$Xmu, k = k, datatype = "iid")
      class(returnlist2) <- c("tbss", "bss") 
      return(returnlist2)
    }
    
    xtFOBI <- tFOBI(x)
    
    x <- xtFOBI$S
    rotat <- tJADERotate(x, k, maxiter = maxiter, eps = eps)
    x <- rotat$x
    
    W <- vector("list", r)
    for(i in 1:r){
      if(k[i] == 0){
        W[[i]] <- diag(dim(x)[i])
        x <- tensorTransform(x, solve(xtFOBI$W[[i]]), i)
      }
      else{
        W[[i]] <- rotat$U[[i]]%*%xtFOBI$W[[i]]
      }
    }
    
    returnlist <- list(S = x, W = W, Xmu = xtFOBI$Xmu, k = k, datatype = "iid")
    
    class(returnlist) <- c("tbss", "bss") 
    
    return(returnlist)
  }


