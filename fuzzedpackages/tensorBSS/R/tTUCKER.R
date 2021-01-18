
tTUCKER <- function(x, ranks, maxiter=1000, eps=1e-06)
{
  r <- length(dim(x)) - 1
  xmu <- apply(x, 1:r, mean)
  xc <- tensorCentering(x)
  init <- tPCA(xc, d=ranks)
  
  
  Us1 <- Us0 <- init$U
  Diff <- Inf
  iter <- 1
  
  while (Diff > eps )
  {
    for (i in 1:r)
    {
      y <- tensorTransform2(xc, Us0, -i, transpose=TRUE)
      COVi <- mModeCovariance(y, i, center = FALSE)
      EV.COVi <- eigen(COVi, symmetric = TRUE)
      Us1[[i]] <- EV.COVi$vectors[, 1:ranks[i], drop = FALSE]
    }
    iter <- iter + 1
    Diff <- sum(unlist(Map("-", lapply(Us1,abs), lapply(Us0,abs)))^2)
    Us0 <- Us1
    if (iter > maxiter) stop("maxiter reached without convergence")
    #print(paste0("iter:",iter," Diff:",Diff))
  }
  y <- tensorTransform2(xc, Us0, 1:r, transpose=TRUE)
  
  
  norm2xc <- sum(xc^2)
  norm2rxc <- sum(tensorTransform2(y, Us1)^2)
  ratio <- norm2rxc/norm2xc
  
  mEV <- vector("list",r)
  names(mEV) <- paste0("Mode",1:r)
  
  for (i in 1:r)
  {
    data.i <- tensorTransform2(xc, Us1, -i, transpose=TRUE)
    mEV[[i]] <- eigen(mModeCovariance(data.i, i, center = TRUE), symmetric = TRUE, only.values = TRUE)$values
  }
  
  res <- list(S=y, U=Us1, Xmu=xmu, tPCA = init, 
              norm2xc = norm2xc,
              norm2rxc = norm2rxc, 
              norm2ratio = ratio,
              mEV = mEV
              )
  res
  
}