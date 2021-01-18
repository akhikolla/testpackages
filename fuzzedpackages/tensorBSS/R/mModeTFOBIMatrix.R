mModeTFOBIMatrix <-
function(x, m, norm){
  xm <- mFlatten(x, m)
  if(norm == FALSE){
    return(mFOBIMatrix(xm))
  }
  else{
    return(mFOBIMatrixNorm(xm))
  }
}
