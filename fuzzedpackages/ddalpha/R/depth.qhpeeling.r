################################################################################
# File:             depth.qhpeeling.r
# Created by:       Pavlo Mozharovskyi
# First published:  08.03.2018
# Last revised:     08.03.2018
# 
# Computation of the convex hull peeling depth.
################################################################################

.qhpeeling_validate <- function(ddalpha, ...) {
  return (list())
}

.qhpeeling_learn <- function(ddalpha) {
  # Calculate depths for each pattern
  for (i in 1:ddalpha$numPatterns){
    ddalpha$patterns[[i]]$depths = 
      .qhpeeling_depths(ddalpha, ddalpha$patterns[[i]]$points)
  }
  # Return the updated structure
  return (ddalpha)
}

.qhpeeling_depths <- function(ddalpha, objects){
  depths <- NULL
  objects = data.matrix(objects)
  for (j in 1:ddalpha$numPatterns){
    pattern <- as.matrix(ddalpha$patterns[[j]]$points)
    
    ds <- rep(0, nrow(objects))
    tmpData <- pattern
    for (i in 1:nrow(pattern)){
      if (length(tmpData) < ncol(pattern) * (ncol(pattern) + 1) + 0.5){
        break
      }
      ds <- ds + as.vector(is.in.convex(objects, tmpData, nrow(tmpData)))
      tmpData <- tmpData[-unique(as.vector(convhulln(tmpData))),]
    }
    
    ds = ds / nrow(pattern)
    
    depths <- cbind(depths, ds)
  }
  return (depths)
}

depth.qhpeeling <- function(x, data){
  if (!is.matrix(x) 
      && is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.numeric(x)){
    stop("Argument \"x\" should be numeric")
  }
  
  if (ncol(x) != ncol(data)){
    stop("Dimensions of the arguments \"x\" and \"data\" should coincide")
  }
  if (ncol(data) + 1 > nrow(data)){ #?
    stop("Too few data points")
  }
  if (is.data.frame(x)){
    x <- as.matrix(x)
  }
  if (is.data.frame(data)){
    data <- as.matrix(data)
  }

  depths <- rep(0, nrow(x))
  tmpData <- data
  for (i in 1:nrow(data)){
    if (length(tmpData) < ncol(data) * (ncol(data) + 1) + 0.5){
      break
    }
    depths <- depths + as.vector(is.in.convex(x, tmpData, nrow(tmpData)))
    tmpData <- tmpData[-unique(as.vector(convhulln(tmpData))),]
  }
  
  return (depths / nrow(data))
}
