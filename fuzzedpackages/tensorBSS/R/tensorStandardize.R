tensorStandardize <-
function(x, location = NULL, scatter = NULL){
  r <- length(dim(x)) - 1
  
  x <- tensorCentering(x, location)
  
  if(is.null(scatter)){
    scatter <- vector("list", r)
    for(m in 1:r){
      scatter[[m]] <- mModeCovariance(x, m, center = FALSE)
    }
  }
  
  invCovStack <- vector("list", r)
  for(m in 1:r){
    invCovStack[[m]] <- symmetricPower(scatter[[m]], -0.5)
  }
  
  for(m in 1:r){
    x <- tensorTransform(x, invCovStack[[m]], m)
  }
  
  attr(x, "location") <- location
  attr(x, "scatter") <- scatter
  
  return(list(x = x, S = invCovStack))
}
