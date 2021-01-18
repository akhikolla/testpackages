################################################################################
# File:             depth.simplicialVolume.r
# Created by:       Oleksii Pokotylo, Pavlo Mozharovskyi
# First published:  15.06.2015
# Last revised:     20.02.2019
# 
# Computation of the simplicial volume data depth.
################################################################################

.longtoint <- function(k){
  limit = 2000000000
  k1 = as.integer(k/limit)
  k2 = k - k1*limit
  return(c(k1, k2))
}

depth.simplicialVolume <- function(x, data, exact = F, k = 0.05, mah.estimate = "moment", mah.parMcd = 0.75, seed = 0){
  if (seed!=0) set.seed(seed)
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
    stop("To few data points")
  }
  if (!exact) if (k <= 0) stop("k must be positive")
      else if (k < 1) k = choose(nrow(data), ncol(data))*k
  if(toupper(mah.estimate) == "NONE"){
    useCov <- 0
    covEst <- diag(ncol(data))
  } else if(toupper(mah.estimate) == "MOMENT"){
    useCov <- 1
    covEst <- cov(data)
  } else if(toupper(mah.estimate) == "MCD"){
    useCov <- 2
    covEst <- covMcd(data, mah.parMcd)$cov
  } else {stop("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")}
      
  points <- as.vector(t(data))
  objects <- as.vector(t(x))
  ds <- .C("OjaDepth", 
           as.double(points), 
           as.double(objects), 
           as.integer(nrow(data)), 
           as.integer(nrow(x)), 
           as.integer(ncol(data)), 
           as.integer(seed),
           as.integer(exact),
           as.integer(.longtoint(k)),
           as.integer(useCov),
           as.double(as.vector(t(covEst))),
           depths=double(nrow(x)))$depths
  
  return (ds)
}



depth.space.simplicialVolume <- function(data, cardinalities, exact = F, k = 0.05, mah.estimate = "moment", mah.parMcd = 0.75, seed = 0){
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.vector(cardinalities, mode = "numeric") 
      || is.na(min(cardinalities)) 
      || sum(.is.wholenumber(cardinalities)) != length(cardinalities) 
      || min(cardinalities) <= 0 
      || sum(cardinalities) != nrow(data)){
    stop("Argument \"cardinalities\" should be a vector of cardinalities of the classes in \"data\" ")
  }
  if (sum(cardinalities < ncol(data) + 1) != 0){
    stop("Not in all classes sufficiently enough objetcs")
  }
  if (toupper(mah.estimate) == "NONE"){
    useCov <- 0
  } else if (toupper(mah.estimate) == "MOMENT"){
    useCov <- 1
  } else if (toupper(mah.estimate) == "MCD"){
    useCov <- 2
  } else {stop("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")}
  
  depth.space <- NULL
  for (i in 1:length(cardinalities)){
    pattern <- data[(1 + sum(cardinalities[0:(i - 1)])):sum(cardinalities[1:i]),]
    pattern.depths <- depth.simplicialVolume(data, pattern, exact, k, mah.estimate, mah.parMcd, seed)
    depth.space <- cbind(depth.space, pattern.depths, deparse.level = 0)
  }
  
  return (depth.space)
}
