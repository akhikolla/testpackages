generate_Gaussian <- function(n,R,p,l,s,sparsityA,sparsityobserved,sparsitylatent,lwb,upb,seed){
  #require(MASS)
  #require(glasso)

  set.seed(seed)
  # Create inverse covariance

  # Generate an Erdos Renyi type network with positive and negative entries
  sparse <- rbinom((p+l)*(p+l),1,1-sparsitylatent)*runif((p+l)*(p+l),lwb,upb)
  #sparse <- rbinom(p*p,1,1-sparsity)*sample(c(-1,1),p*p,replace=TRUE)*runif(p*p,lwb,upb)
  Theta <- matrix(data=sparse,nrow=p+l,ncol=p+l)
  Theta[lower.tri(Theta,diag=FALSE)] <- 0
  Theta <- Theta+t(Theta)
  diag(Theta) <- 0

  invSigma <- Theta

  # Generate an Erdos Renyi type network with positive and negative entries
  sparse <- rbinom(p*p,1,1-sparsityobserved)*runif(p*p,lwb,upb)
  #sparse <- rbinom(p*p,1,1-sparsity)*sample(c(-1,1),p*p,replace=TRUE)*runif(p*p,lwb,upb)
  Theta <- matrix(data=sparse,nrow=p,ncol=p)
  Theta[lower.tri(Theta,diag=FALSE)] <- 0
  Theta <- Theta+t(Theta)
  diag(Theta) <- 0

  invObs <- Theta

  invSigma[1:p,1:p] <- invObs
  ee <- min(eigen(invSigma,only.values=T)$values)
  diag(invSigma) <- diag(invSigma)+ifelse(ee < 0, -ee + 0.2, 0.2)
  a <- matrix(0,nrow=p+l,ncol=p+l)
  diag(a) <- diag(invSigma)
  invSigma <- sqrt(solve(a))%*%invSigma%*%sqrt(solve(a))
  Sigma <- solve(invSigma)

  # Generate latent variable for the n observations
  latent <- mvrnorm(n*s,rep(0,p+l),Sigma)
  # latentnew<-latent[rep(1:nrow(latent),c(m,n,p,...q),]
  t <- R/s
  latentnew <- latent[rep(1:nrow(latent),each=t),]


  # conditional distribution of observed given latent
  partialSigma <- Sigma[1:p,1:p]-Sigma[1:p,(p+1):(p+l)]%*%solve(Sigma[(p+1):(p+l),(p+1):(p+l)])%*%Sigma[(p+1):(p+l),1:p]

  # generate A
  #sparse <- rbinom(p*p,1,1-sparsityA)*runif(p*p,0.4,0.4)
  sparse <- rbinom(p*p,1,1-sparsityA)*sample(c(-1,1),p*p,replace=TRUE)*runif(p*p,0.3,0.3)
  A <- matrix(data=sparse,nrow=p,ncol=p)
  diag(A) <- 0.3
  #diag_a <- sample(c(-1,1),p,replace=TRUE)*runif(p,0.45,0.45)
  #diag(A) <- diag_a

  dat <- NULL
  for(i in 1:n){
    partialmu <- Sigma[1:p,(p+1):(p+l)]%*%solve(Sigma[(p+1):(p+l),(p+1):(p+l)])%*%latentnew[1+R*(i-1),(p+1):(p+l)]
    Y <- matrix(mvrnorm(1,partialmu,partialSigma), nrow = 1, ncol = p)
    for(j in 1:(R-1)){
      partialmu <- Sigma[1:p,(p+1):(p+l)]%*%solve(Sigma[(p+1):(p+l),(p+1):(p+l)])%*%latentnew[1+R*(i-1)+j,(p+1):(p+l)]
      mu <-  A%*%Y[j,] + partialmu
      Y <- rbind(Y,mvrnorm(1,mu,partialSigma))
    }
    dat <- rbind(dat,Y)
  }

  X <- array(data=NA,c(R,p,n))
  for(i in 1:n){
    X[,,i] <- dat[(1+(i-1)*R):(i*R),]
  }

  # Transform the data X into a list
  Xtemp <- aperm(X,c(3,1,2))
  X <- vector('list',n)
  for(i in 1:n){
    X[[i]] <- Xtemp[i,,]
  }

  # true graph
  inv_PartialSigma <- solve(partialSigma)
  inv_PartialSigma <- ifelse(abs(inv_PartialSigma) < 1e-7, 0, inv_PartialSigma)
  diag(inv_PartialSigma) <- 0

  return(list(X = X, truegraph = inv_PartialSigma))
}
