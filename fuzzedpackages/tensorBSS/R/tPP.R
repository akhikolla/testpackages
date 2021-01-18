tPP <-
  function(x, nl = "pow3", eps = 1e-6, maxiter = 100){
    
    if(length(dim(x)) == 2){
      stop("Does not currently work for vector-valued data. Try the package ICtest and the function NGPP.")
    }
    
    # Parse non-linearities to integers
    nl_int <- match(nl, c("pow3", "skew", "tanh"))
    
    if(nl_int == 3){
      stop("The non-linearity tanh is not currently supported.")
    }
    
    # Needed for TFOBI-rotations
    norm <- rep(FALSE, length(dim(x)) - 1)
    
    xmu <- apply(x, 1:(length(dim(x)) - 1), mean)
    
    stand <- tensorStandardize(x)
    xst <- stand$x
    rotat <- tFOBIRotate(xst, norm)
    
    r <- length(stand$S)
    W <- vector("list", r)
    s <- xst
    iter <- rep(0, r)
    for(i in 1:r){
      U0 <- t(rotat$U[[i]])
      # p <- dim(x)[i]
      # U0 <- rorth(p)
      rotat_2 <- tPPRotate(U0, mFlatten(xst, i), nl_int, eps, maxiter)
      U0 <- rotat_2[[1]]
      iter[i] <- rotat_2[[2]]
      
      W[[i]] <- t(U0)%*%stand$S[[i]]
      s <- tensorTransform(s, t(U0), i) 
    }
    
    returnlist <- list(S = s, W = W, iter = iter, Xmu = xmu, datatype = "iid")
    
    class(returnlist) <- c("tbss", "bss") 
    
    return(returnlist)
  }
