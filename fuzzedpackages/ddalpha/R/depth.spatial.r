depth.spatial.local <- function(x, data, mah.estimate = "moment", mah.parMcd = 0.75, kernel.bandwidth = 1){
  if (!is.matrix(x) 
      && is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  mean <- colMeans(data)
  if(mah.estimate == "none"){
    lambda = diag(ncol(data))
  } else {
    if(mah.estimate == "moment"){
      cov <- cov(data)
    } else if(mah.estimate == "MCD"){
      cov <- covMcd(data, mah.parMcd)$cov
    } else {stop("Wrong parameter 'mah.estimate'")}
    cov.eig <- eigen(cov)
    B <- cov.eig$vectors %*% diag(sqrt(cov.eig$values))
    lambda <- solve(B)
  }
  
  if (kernel.bandwidth>1){
    k <- function(x){
      return(sqrt(2*pi)^(-ncol(x))*exp(-rowSums(x^2)/2/kernel.bandwidth))
    }
  }else{  #1/kernel.bandwidth^d
    k <- function(x){
      return((sqrt(2*pi)*kernel.bandwidth)^(-ncol(x))*exp(-rowSums(x^2)/2/kernel.bandwidth))
    }
  }  

  depths <- apply(x, 1, function(x){  
            n = nrow(data)
            t1 = t(lambda %*% (x - t(data)))
            t1 <- t1[which(rowSums(t1) != 0),]
            
            return (sum(k(t1))/n - sqrt(sum((colSums( k(t1)*( t1/sqrt(rowSums(t1^2)) ) )/n)^2)))
          }  
        )
  return (depths)
}

depth.spatial <- function(x, data, mah.estimate = "moment", mah.parMcd = 0.75){
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
  mean <- colMeans(data)
  if(mah.estimate == "none"){
    lambda = diag(ncol(data))
  } else {
    if(mah.estimate == "moment"){
      cov <- cov(data)
    } else if(mah.estimate == "MCD"){
      cov <- covMcd(data, mah.parMcd)$cov
    } else {stop("Wrong parameter 'mah.estimate'")}
    if(sum(is.na(cov)) == 0){
      cov.eig <- eigen(cov)
      B <- cov.eig$vectors %*% diag(sqrt(cov.eig$values))
      lambda <- solve(B)
    } else{
      lambda = diag(ncol(data))
    }
  }
  
  depths <- rep(-1, nrow(x))
  for (i in 1:nrow(x)){
    tmp1 <- t(lambda %*% (x[i,] - t(data)))
    tmp1 <- tmp1[which(rowSums(tmp1) != 0),]
    tmp2 <- 1/sqrt(rowSums(tmp1^2))
    depths[i] <- 1 - sqrt(sum((colSums(tmp2*tmp1)/nrow(data))^2))
  }
  return (depths)
}

depth.space.spatial <- function(data, cardinalities, mah.estimate = "moment", mah.parMcd = 0.75){
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
  
  depth.space <- NULL
  for (i in 1:length(cardinalities)){
    pattern <- data[(1 + sum(cardinalities[0:(i - 1)])):sum(cardinalities[1:i]),]
    pattern.depths <- depth.spatial (data, pattern, mah.estimate, mah.parMcd)
    depth.space <- cbind(depth.space, pattern.depths, deparse.level = 0)
  }
  
  return (depth.space)
}