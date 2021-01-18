W.mat <- function(pi) {
  nind <- nrow(pi)
  ngpe <- ncol(pi)
  
  M <- vector("list", ngpe)
  for(i in 1:ngpe) M[[i]] <- do.call(rbind, lapply(1:ngpe, function(j) diag(-pi[,i]*pi[,j])))
  W <- do.call(cbind, M)
  
  diag(W) <- as.vector(pi*(1-pi))
  
  W
}

P.mat <- function(pi, X) {
  ngpe <- ncol(pi)
  W <- W.mat(pi)
  XX <- .Call("block_diag", PACKAGE = "Ravages", rep(list(X), ngpe) )
  WX <- W %*% XX
  XWX <- t(XX) %*% WX
  W - WX %*% solve( XWX, t(WX) )
}

Gamma.mat <- function(nind, ngpe) {
  Id.bygpe <- diag(1, ncol=nind, nrow=nind)
  Gamma.prov <- .Call("block_diag", PACKAGE = "Ravages", rep(list(Id.bygpe), ngpe-1))
  mId.gpe1 <- do.call(cbind, rep(list(diag(-1, nrow=nind)), ngpe-1))
  rbind(mId.gpe1, Gamma.prov)
}

P.mat2 <- function( pi, X ) {
  P <- P.mat( pi[,-1, drop = FALSE], X)
  Gamma <- Gamma.mat(nrow(pi), ncol(pi))
  Gamma %*% P %*% t(Gamma)
}

#Not used: new C++ function
block.diag <- function(L) {
  if(length(L) == 1)
    return(L[[1]])
  D <- block.diag( L[-1] )
  L1 <- L[[1]]
  ZerosNE <- matrix(0, nrow = nrow(L1), ncol = ncol(D))
  ZerosSW <- matrix(0, nrow = nrow(D), ncol(L1))
  rbind( cbind(L1, ZerosNE),   cbind(ZerosSW, D) )
}
