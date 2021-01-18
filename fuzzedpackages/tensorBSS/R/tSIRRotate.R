tSIRRotate <-
  function(x, y, h, ...){
    r <- length(dim(x)) - 1
    
    rotateStack <- vector("list", r)
    for(m in 1:r){
      rotateStack[[m]] <- t(eigenVectors(mModeTSIRMatrix(x, m, y, h, ...)))
    }
    
    
    for(m in 1:r){
      x <- tensorTransform(x, rotateStack[[m]], m)
    }
    
    return(list(x = x, U = rotateStack))
  }
