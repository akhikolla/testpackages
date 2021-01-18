semilatent <- function(data, n, R, p, lambda, distribution = "Gaussian", rule = "AND"){
  if (class(data)[1] == "matrix") {
    dat <- array(data=NA,c(R,p,n))
    for(i in 1:n){
      dat[,,i] <- data[(1+(i-1)*R):(i*R),]
    }
    dat_temp <- aperm(dat,c(3,1,2))
    dat <- vector('list',n)
    for(i in 1:n){
      dat[[i]] <- dat_temp[i,,]
    }
  }

  if (class(data)[1] == "list") {
    dat <- data
  }

  if (class(data)[1] == "array") {
    dat <- list(data[,,1])
    for (i in 2:n) {
      dat[[i]] <- data[,,i]
    }
  }

  if (class(data)[1] == "data.frame") {
    data <- as.matrix(data)
    dat <- array(data=NA,c(R,p,n))
    for(i in 1:n){
      dat[,,i] <- data[(1+(i-1)*R):(i*R),]
    }
    dat_temp <- aperm(dat,c(3,1,2))
    dat <- vector('list',n)
    for(i in 1:n){
      dat[[i]] <- dat_temp[i,,]
    }
  }

  newdata <- transformDataReplicates(dat,n,p,rep(R,n))

  if (distribution == "Gaussian") {
    omega <- semigraph_GGM(newdata,n,p,lambda)$beta
  }

  if (distribution == "Ising") {
    omega <- semigraph_Ising(newdata,n,p,lambda)$beta
  }
  
  res <- omega
  # AND rule
  if (rule == "AND") {
    res <- res!=0
    res <- res*1
    res <- res + t(res)
    res[which(res==1)] <- 0
  }
  # OR rule
  if (rule == "OR") {
    res <- (res+t(res))/2
  }
  
  rr <- res!=0
  rr <- rr*1
  return(list(omega=omega, theta = rr, penalty = lambda))
}
