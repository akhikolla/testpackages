tFOBIRotate <-
function(x, norm){
  r <- length(dim(x)) - 1
  
  rotateStack <- vector("list", r)
  for(m in 1:r){
    rotateStack[[m]] <- t(eigenVectors(mModeTFOBIMatrix(x, m, norm[m])))
  }
  
  
  for(m in 1:r){
    x <- tensorTransform(x, rotateStack[[m]], m)
  }
  
  return(list(x = x, U = rotateStack))
}
