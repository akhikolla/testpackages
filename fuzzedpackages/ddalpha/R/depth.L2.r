################################################################################
# File:             depth.L2.r
# Created by:       Pavlo Mozharovskyi
# First published:  08.03.2018
# Last revised:     08.03.2018
# 
# Computation of the L2-depth.
################################################################################

.L2_validate <- function(ddalpha, mah.estimate = "moment", mah.parMcd = 0.75, ...) {
  # only validate and stop if anything is wrong
    if( !(toupper(mah.estimate) %in% c("NONE", "MOMENT", "MCD"))){
      stop("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")
    }

  return (list(mah.estimate = mah.estimate, mah.parMcd = mah.parMcd))
}

.L2_learn <- function(ddalpha) {
  #1. Calculate statistics based on data
  
  if(toupper(ddalpha$mah.estimate) == "NONE"){
    for (i in 1:ddalpha$numPatterns){
      ddalpha$patterns[[i]]$sigma <- diag(ddalpha$dimension)
    }
  } else {
    for (i in 1:ddalpha$numPatterns){
      if(toupper(ddalpha$mah.estimate) == "MOMENT"){
        cov <- cov(as.matrix(ddalpha$patterns[[i]]$points))
      } else if(toupper(ddalpha$mah.estimate) == "MCD"){
        cov <- covMcd(as.matrix(ddalpha$patterns[[i]]$points), ddalpha$mah.parMcd)$cov
      } else {stop("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")}
      if(sum(is.na(cov)) == 0){
        ddalpha$patterns[[i]]$sigma <- solve(cov)
      } else{
        ddalpha$patterns[[i]]$sigma <- diag(ddalpha$dimension)
        warning("Covariance estimate not found for pattern ", ddalpha$patterns[[i]]$name, ", no affine-invariance-adjustment")
      }
    }
  }
  
  #2. Calculate depths for each pattern
  for (i in 1:ddalpha$numPatterns){
    ddalpha$patterns[[i]]$depths = .L2_depths(ddalpha, ddalpha$patterns[[i]]$points)
  }
  
  return (ddalpha)
}

.L2_depths <- function(ddalpha, objects){
  depths <- NULL
  objects = data.matrix(objects)
  for (j in 1:ddalpha$numPatterns){
    pattern <- as.matrix(ddalpha$patterns[[j]]$points)

    sigma = ddalpha$patterns[[j]]$sigma
    
    ds <- rep(-1, nrow(objects))
    for (i in 1:nrow(objects)){
      tmp1 <- t(objects[i,] - t(pattern))
      tmp2 <- tmp1 %*% sigma
      ds[i] <- 1/(1 + mean(sqrt(rowSums(tmp2 * tmp1))))
    }
    
    depths <- cbind(depths, ds)
  }
  return (depths)
}

depth.L2 <- function(x, data, mah.estimate = "moment", mah.parMcd = 0.75){
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if(is.data.frame(data))
    data = data.matrix(data)
  if (!is.matrix(x)){
    if(is.vector(x))
      x <- matrix(x, nrow=1)
    if(is.data.frame(x))
      x = data.matrix(x)
  }
  if(toupper(mah.estimate) == "NONE"){
    sigma = diag(ncol(data))
  } else {
    if(toupper(mah.estimate) == "MOMENT"){
      cov <- cov(data)
    } else if(toupper(mah.estimate) == "MCD"){
      cov <- covMcd(data, mah.parMcd)$cov
    } else {stop("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")}
    if(sum(is.na(cov)) == 0){
      sigma <- solve(cov)
    } else{
      sigma = diag(ncol(data))
      warning("Covariance estimate not found, no affine-invariance-adjustment")
    }
  }

  depths <- rep(-1, nrow(x))
  for (i in 1:nrow(x)){
    tmp1 <- t(x[i,] - t(data))
    tmp2 <- tmp1 %*% sigma
    depths[i] <- 1/(1 + mean(sqrt(rowSums(tmp2 * tmp1))))
  }
  return (depths)
}
