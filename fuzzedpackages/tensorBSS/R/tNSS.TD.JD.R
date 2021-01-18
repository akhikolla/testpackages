tNSS.TD.JD <- function(x, K = 12, lags = 0:12, n.cuts = NULL, eps = 1e-06, maxiter = 100, ...)
{
  
  dim_x <- dim(x)
  r <- length(dim_x) - 1
  n <- dim_x[r + 1]
  
  if(length(dim_x) == 2){
    Xmu <- rowMeans(x)
    returnlist <- NSS.TD.JD(t(x), K = K, Tau = lags, n.cuts = n.cuts, eps = eps, maxiter = maxiter, ...)
    returnlist$S <- t(returnlist$S)
    returnlist2 <- list(S = returnlist$S,
                        W = returnlist$W,
                        K = returnlist$K,
                        lags = returnlist$k,
                        n.cuts = returnlist$n.cut,
                        Xmu = Xmu,
                        datatype = "ts")
    class(returnlist2) <- c("tbss", "bss") 
    return(returnlist2)
  }
  
  if(is.null(n.cuts))
  {
    if(K == 1)
    {
      slices <- rep(1, n)
      n.cuts <- c(0, n)
    } else {
      slices <- as.numeric(cut(1:n, breaks = K, labels = 1:K))
      n.cuts <- c(0, which(slices[-n] - slices[-1] == -1), n)
    }
  } else {
    K <- length(n.cuts) - 1
    slices <- as.numeric(cut(1:n, breaks = n.cuts, labels = 1:K))
  }
  
  
  # Save the mean
  xmu <- apply(x, 1:r, mean)
  
  # Standardize the data
  res_stand <- tensorStandardize(x)
  x_stand <- res_stand$x
  
  U_list <- vector("list", r)
  
  # Iterate over modes
  for(m in 1:r)
  {
    current_dim <- dim_x[m]
    matrix_array <- array(0, dim=c(current_dim, current_dim, K*length(lags)))
    
    # Iterate over slices
    for(h in 1:K)
    {
      
      # extractor <- diag(n)[slices == h, ]
      # x_slice <- aperm(tensorTransform(aperm(x_stand), extractor, 1))
      
      # Select the obs. belonging to the current slice
      x_slice <- arraySelectLast(x_stand, (slices == h))
      
      # Iterate over lags
      for(tau in lags){
        
        # Compute the matrix of interest (with re-centering within intervals)
        mMAC <- mModeAutoCovariance(x_slice, m, tau, center = TRUE)
        matrix_array[, , (h - 1)*length(lags) + match(tau, lags)] <- 0.5*(mMAC + t(mMAC))
      }
    }
    
    U_list[[m]] <- frjd(matrix_array, eps = eps, maxiter = maxiter, ...)$V
    
  }
  
  # Calculate the components
  S <- x_stand
  
  for(m in 1:r)
  {
    S <- tensorTransform(S, t(U_list[[m]]), m)
  }
  
  # Calculate the unmixing matrices
  
  W <- vector("list", r)
  for(m in 1:r)
  {
    W[[m]] <- crossprod(U_list[[m]], res_stand$S[[m]])
  }
  
  returnlist <- list(S = S,
                     W = W,
                     K = K,
                     lags = lags,
                     n.cuts = n.cuts,
                     Xmu = xmu,
                     datatype = "ts")
  class(returnlist) <- c("tbss", "bss") 
  
  return(returnlist)
}
