tNSS.SD <- function(x, n.cuts = NULL)
{
  
  dim_x <- dim(x)
  r <- length(dim_x) - 1
  n <- dim_x[r + 1]
  
  if(length(dim_x) == 2){
    Xmu <- rowMeans(x)
    returnlist <- NSS.SD(t(x), n.cut = n.cuts)
    returnlist$S <- t(returnlist$S)
    returnlist2 <- list(S = returnlist$S,
                        W = returnlist$W,
                        EV = returnlist$EV,
                        n.cuts = returnlist$n.cut,
                        Xmu = Xmu,
                        datatype = "ts")
    class(returnlist2) <- c("tbss", "bss") 
    return(returnlist2)
  }
  
  if(any(is.null(n.cuts))){
    n.cuts <- c(0, floor((n + 1)/2), n)
  }
  
  # Slice the time axis
  slices <- as.numeric(cut(1:n, breaks = n.cuts, labels = 1:2))
  
  # Save the mean
  Xmu <- apply(x, 1:r, mean)
  
  # Center the data
  x_center <- tensorCentering(x)
  
  W_list <- vector("list", r)
  
  # Iterate over modes
  for(m in 1:r)
  {
    
    x_slice_1 <- arraySelectLast(x_center, (slices == 1))
    x_slice_2 <- arraySelectLast(x_center, (slices == 2))
    
    mMAC_1 <- mModeCovariance(x_slice_1, m, center = TRUE)
    symm_mMAC_1 <- 0.5*(mMAC_1 + t(mMAC_1))
    mMAC_2 <- mModeCovariance(x_slice_2, m, center = TRUE)
    symm_mMAC_2 <- 0.5*(mMAC_2 + t(mMAC_2))
    
    eig_mMAC_1 <- eigen(mMAC_1, symmetric = TRUE)
    mMAC_1_root <- eig_mMAC_1$vectors%*%diag(eig_mMAC_1$values^(-1/2))%*%t(eig_mMAC_1$vectors)
    
    eig_final <- eigen(tcrossprod(crossprod(mMAC_1_root, mMAC_2), mMAC_1_root), symmetric = TRUE)
    U <- eig_final$vectors
    
    W_list[[m]] <- crossprod(U, mMAC_1_root)
    
  }
  
  # Calculate the components
  S <- x_center
  
  for(m in 1:r)
  {
    S <- tensorTransform(S, W_list[[m]], m)
  }
  
  returnlist <- list(S = S,
                     W = W_list,
                     EV = eig_final$values,
                     n.cuts = n.cuts,
                     Xmu = Xmu,
                     datatype = "ts")
  class(returnlist) <- c("tbss", "bss") 
  
  return(returnlist)
}