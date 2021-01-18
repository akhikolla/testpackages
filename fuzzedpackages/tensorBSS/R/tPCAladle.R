# tPCAlade functions

List2dataframe <- function(x)  do.call(cbind.data.frame, x)

fi <- function (EVboot, EVdata, r) 
{
  fni <- numeric(r)
  for (ii in 1:r) {
    fni[ii] <- det(crossprod(EVdata[, 1:ii], EVboot[, 1:ii]))
  }
  1 - abs(fni)
}



Boot_tPCA_ladle <- function(x, U,mL){
  r <- length(dim(x)) - 1
  x <- tensorCentering(x)
  FI <- vector("list",r)
  for (m in 1:r) {
    mCov <- mModeCovariance(x, m, center = FALSE)
    mEVD <- eigen(mCov, symmetric = TRUE)
    FI[[m]] <- fi(mEVD$vectors,U[[m]],mL[m])
  }
  FI
}

tPCAladle <- function(x, n.boot = 200, ncomp = NULL)
{
  data.name <- deparse(substitute(x))
  method <- "tPCA"  
  
  r <- length(dim(x)) - 1
  if (is.null(ncomp)) {ncomp <- dim(x)[1:r]; ncomp <- ifelse(ncomp > 10, floor(ncomp /log(ncomp)), ncomp - 1)}
  
  
  xmu <- apply(x, 1:r, mean)
  x <- tensorCentering(x)
  U <- vector("list", r)
  D <- vector("list", r)
  
  
  for (m in 1:r) {
    mCov <- mModeCovariance(x, m, center = FALSE)
    mEig <- eigen(mCov, symmetric = TRUE)
    U[[m]] <- mEig$vectors
    D[[m]] <- mEig$values
  }
  
  fis <- replicate(n.boot, Boot_tPCA_ladle(tensorBoot(x), U, ncomp), simplify=FALSE)
  
  ResMode <- vector("list", r)
  names(ResMode) <- paste0("Mode", 1:r)
  
  for (ii in 1:r) {
    EVii <- D[[ii]] 
    fn0ii <- c(0,rowMeans(sapply(fis, function(x) x[[ii]])))
    fnii <- fn0ii/(1 + sum(fn0ii))
    phinii <- EVii[1:(ncomp[ii] + 1)]/(1 + sum(EVii[1:(ncomp[ii] +  1)]))
    gnii <- fnii + phinii
    est.kii <- which.min(gnii) - 1
    ResMode[[ii]] <- list(mode = paste("Mode",ii),
                          k = est.kii,
                          fn = fnii,
                          phin = phinii,
                          lambda = EVii[1:(ncomp[ii] + 1)],
                          gn = gnii,
                          comp = 0:ncomp[ii]
    )
  }
  
  s <- tensorTransform2(x, U, 1:r, transpose=TRUE)
  
  RES <- list(U=U, D=D, S=s, ResMode = ResMode, xmu =  xmu, data.name = data.name, method = method)
  class(RES) <- "tladle"
  RES
}


