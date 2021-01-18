tPCA <-
function(x, p = NULL, d = NULL){
  
  r <- length(dim(x)) - 1
  
  xmu <- apply(x, 1:r, mean)
  
  x <- tensorCentering(x)
  
  # Compute rotations
  U <- vector("list", r)
  D <- vector("list", r)
  for(m in 1:r){
    mCov <- mModeCovariance(x, m, center = FALSE)
    mEig <- eigen(mCov, symmetric = TRUE)
    U[[m]] <- mEig$vectors
    D[[m]] <- mEig$values
  }
  
  
  d_comp <- dim(x)[1:r]
  p_comp <- NULL
  
  # Reduce dimension, if wanted
  if(!is.null(p)){
    d_comp <- NULL
    for(m in 1:r){
      csm <- cumsum(D[[m]]/sum(D[[m]]))
      d_comp <- c(d_comp, sum(csm < p[m]) + 1)
    }
  }
  else if(!is.null(d)){
    d_comp <- d
  }
  
  for(m in 1:r){
    p_comp <- c(p_comp, cumsum(D[[m]]/sum(D[[m]]))[d_comp[m]])
  }
  
  for(m in 1:r){
    U[[m]] <- U[[m]][, 1:d_comp[m], drop = FALSE]
    # D[[m]] <- D[[m]][1:d_comp[m]]
  }
  
  # Compute the scores
  for(m in 1:r){
    x <- tensorTransform(x, t(U[[m]]), m)  
  }
  
  returnlist <- list(S = x, U = U, D = D, p_comp = p_comp, Xmu = xmu)
  
  return(returnlist)
}