responseFun2 <- function(eta) {
  q <- length(eta)

  eta.help <- matrix(rep(c(0, eta), each = q + 1), ncol = q + 
                       1)
  eta.help[upper.tri(eta.help)] <- 0
  pi <- cumprod(c(1, exp(eta[-q])))/sum(apply(exp(eta.help), 
                                              1, prod))

  pi
}


sim_fun <- function(model, m, I, k, n, gamma, seed = NULL){

  if(!is.null(seed)){
    set.seed(seed)
  }
  
  RSM <- GPCM <- FALSE
  if(model %in% c("GRSM","RSM")){
    RSM <- TRUE
  }
  if(model %in% c("GRSM","GPCM","2PL")){
    GPCM <- TRUE
  }
  
  
  q <- k-1

  X <- c()
  for(i in 1:m){
    if(i%%2 == 1){
      X <- cbind(X, rnorm(n))
    }else{
      X <- cbind(X, rbinom(n,1,0.5))
    }
  }
  X <- scale(X)
  
  if(!RSM){
    delta <- deltaX <- matrix(round(rnorm(q*I,sd=0.5),2),nrow=I)
    # if(model == "GPCM"){
    #   delta <- deltaX <- delta-0.5
    # }
    alpha <- NA
  }else{
    delta <- round(rnorm(I,sd=0.5),2)
    alpha <- c(0,round(rnorm(q-1,sd=0.5),2))
    deltaX <- t(t(matrix(rep(delta,q),nrow=I))+alpha)
  }
  
  
  if(!GPCM){
    sigma <- rep(1,I)
  }else{
    sigma <- seq(0.7,1,length=I)
  }
  
  theta <- rnorm(n)
  
  
  lin_pred<- c()
  probs <- c()
  y <- c()
  
  for(i in 1:n){
    for(ii in 1:I){
      eta <- sigma[ii] * (theta[i] - deltaX[ii,] -sum(gamma[ii,]*X[i,]))
      
      lin_pred <- rbind(lin_pred, eta)
      pi <- responseFun2(eta)
      if(q==1){
        pi <- exp(eta)/(1+exp(eta))
      }
      probs <- rbind(probs,pi)
      pi <- c(pi,1-sum(pi))

      # if(any(pi<0)){browser()}
      y.sample <- which(rmultinom(1,1,pi)==1) 
      y <- c(y,y.sample)
    }
  }
  
  Y <- matrix(y,byrow=TRUE,nrow=n)
  data.sim <- as.data.frame(cbind(Y,X))
  
  return(list(data=data.sim, theta = theta, alpha = alpha, sigma = sigma, 
              delta = delta, gamma = gamma, lin_pred = lin_pred, probs = probs))
}


sim_fun2 <- function(model, m, I, k, n, gamma, seed = NULL){
  
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  RSM <- GPCM <- FALSE
  if(model %in% c("GRSM","RSM")){
    RSM <- TRUE
  }
  if(model %in% c("GRSM","GPCM","2PL")){
    GPCM <- TRUE
  }
  
  
  q <- k-1
  
  X <- c()
  for(i in 1:m){
    X <- cbind(X, rbinom(n,1,0.5))
  }
  
  
  if(!RSM){
    delta <- deltaX <- matrix(round(rnorm(q*I,sd=0.3),2),nrow=I)
    alpha <- NA
  }else{
    delta <- round(rnorm(I,sd=0.3),2)
    alpha <- round(rnorm(q,sd=0.3),2)
    deltaX <- t(t(matrix(rep(delta,q),nrow=I))+alpha)
  }
  
  
  if(!GPCM){
    sigma <- rep(1,I)
  }else{
    sigma <- seq(0.8,1.2,length=I)
  }
  
  theta <- rnorm(n)
  #
  
  lin_pred<- c()
  probs <- c()
  y <- c()
  
  for(i in 1:n){
    for(ii in 1:I){
      eta <- sigma[ii] * (theta[i] - deltaX[ii,] -sum(gamma[ii,]*X[i,]))
      
      lin_pred <- rbind(lin_pred, eta)
      pi <- responseFun2(eta)
      if(q==1){
        pi <- exp(eta)/(1+exp(eta))
      }
      probs <- rbind(probs,pi)
      pi <- c(pi,1-sum(pi))
      
      # if(any(pi<0)){browser()}
      y.sample <- which(rmultinom(1,1,pi)==1) 
      y <- c(y,y.sample)
    }
  }
  
  Y <- matrix(y,byrow=TRUE,nrow=n)
  data.sim <- as.data.frame(cbind(Y,X))
  
  return(list(data=data.sim, theta = theta, alpha = alpha, sigma = sigma, 
              delta = delta, gamma = gamma, lin_pred = lin_pred, probs = probs))
}


sim_fun3 <- function(model, m, I, k, n, gamma, seed = NULL){
  
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  RSM <- GPCM <- FALSE
  if(model %in% c("GRSM","RSM")){
    RSM <- TRUE
  }
  if(model %in% c("GRSM","GPCM","2PL")){
    GPCM <- TRUE
  }
  
  
  q <- k-1
  
  for(i in 1:m){
    if(i%%3 == 1){
      X <- data.frame(V1= rnorm(n))
    }
    if(i%%3 == 2){
      X$V2 <- factor(rbinom(n,1,0.5))
    }
    if(i%%3 == 0){
      X$V3 <- factor(sample(1:4,n,replace=TRUE))
    }
  }
  
  
  X2 <- model.matrix(~V1+V2+V3,data=X)[,-1]
  X2 <- scale(X2)
  
  if(!RSM){
    delta <- deltaX <- matrix(round(rnorm(q*I,sd=0.5,mean=-0.5),2),nrow=I)
    # if(model == "GPCM"){
    #   delta <- deltaX <- delta-0.5
    # }
    alpha <- NA
  }else{
    delta <- round(rnorm(I,sd=0.5),2)
    alpha <- c(0,round(rnorm(q-1,sd=0.5),2))
    deltaX <- t(t(matrix(rep(delta,q),nrow=I))+alpha)
  }
  
  
  if(!GPCM){
    sigma <- rep(1,I)
  }else{
    sigma <- seq(0.7,1,length=I)
  }
  
  theta <- rnorm(n)
  
  lin_pred<- c()
  probs <- c()
  y <- c()
  
  for(i in 1:n){
    for(ii in 1:I){
      eta <- sigma[ii] * (theta[i] - deltaX[ii,] -sum(gamma[ii,]*X2[i,]))
      
      lin_pred <- rbind(lin_pred, eta)
      pi <- responseFun2(eta)
      if(q==1){
        pi <- exp(eta)/(1+exp(eta))
      }
      probs <- rbind(probs,pi)
      pi <- c(pi,1-sum(pi))
      
      # if(any(pi<0)){browser()}
      y.sample <- which(rmultinom(1,1,pi)==1) 
      y <- c(y,y.sample)
    }
  }
  
  Y <- matrix(y,byrow=TRUE,nrow=n)
  data.sim <- as.data.frame(cbind(Y,X))
  names(data.sim)[1:I] <- paste0("Item",1:I)
  
  return(list(data=data.sim, theta = theta, alpha = alpha, sigma = sigma, 
              delta = delta, gamma = gamma, lin_pred = lin_pred, probs = probs))
}



sim_cor <- function(model, m, I, k, n, gamma, sigma, seed = NULL){
  
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  RSM <- GPCM <- FALSE
  if(model %in% c("GRSM","RSM")){
    RSM <- TRUE
  }
  if(model %in% c("GRSM","GPCM","2PL")){
    GPCM <- TRUE
  }
  
  
  q <- k-1
  
  mat1 <- rmvnorm(n, sigma = sigma)
  
  set.seed(1860)
  
  X <- c()

  for(i in 1:m){
    if(i%%2 == 1){
      X <- cbind(X, mat1[,i])
    }else{
      X <- cbind(X, mat1[,i]>0)
    }
  }
  X <- scale(X)
  
  if(!RSM){
    delta <- deltaX <- matrix(round(rnorm(q*I,sd=0.5),2),nrow=I)
    # if(model == "GPCM"){
    #   delta <- deltaX <- delta-0.5
    # }
    alpha <- NA
  }else{
    delta <- round(rnorm(I,sd=0.5),2)
    alpha <- c(0,round(rnorm(q-1,sd=0.5),2))
    deltaX <- t(t(matrix(rep(delta,q),nrow=I))+alpha)
  }
  
  
  if(!GPCM){
    sigma <- rep(1,I)
  }else{
    sigma <- seq(0.7,1,length=I)
  }
  
  theta <- rnorm(n)
  
  
  lin_pred<- c()
  probs <- c()
  y <- c()
  
  for(i in 1:n){
    for(ii in 1:I){
      eta <- sigma[ii] * (theta[i] - deltaX[ii,] -sum(gamma[ii,]*X[i,]))
      
      lin_pred <- rbind(lin_pred, eta)
      pi <- responseFun2(eta)
      if(q==1){
        pi <- exp(eta)/(1+exp(eta))
      }
      probs <- rbind(probs,pi)
      pi <- c(pi,1-sum(pi))
      
      # if(any(pi<0)){browser()}
      y.sample <- which(rmultinom(1,1,pi)==1) 
      y <- c(y,y.sample)
    }
  }
  
  Y <- matrix(y,byrow=TRUE,nrow=n)
  data.sim <- as.data.frame(cbind(Y,X))
  
  return(list(data=data.sim, theta = theta, alpha = alpha, sigma = sigma, 
              delta = delta, gamma = gamma, lin_pred = lin_pred, probs = probs))
}
