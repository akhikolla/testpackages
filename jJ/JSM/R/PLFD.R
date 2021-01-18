
#================== Profile Likelihood Method with Forward Difference for Model I & II ==================#

PLFD <- function (model, theta, cvals, n, ncx, ncz, ncw, varNames, p, Z.st, Y.st, X.st, b, Ztime, nk, Wtime, Xtime, Wtime2, Xtime2, GQ, rho, d, wGQ, Ztime2.st, Indcs) {

  delta = cvals$delta;
  tol = min(cvals$tol.P, cvals$delta)/100;
  iter = cvals$max.iter
  alpha.name = varNames$alpha.name;
  phi.names = varNames$phi.names
  beta.names = varNames$beta.names

  Index = Indcs$Index  
  Index0 = Indcs$Index0
  Index1 = Indcs$Index1
  Index2 = Indcs$Index2


  pl <- LHGeneric(theta, n, Z.st, Y.st, X.st, b, Ztime, nk, Index0, model, Wtime, Xtime, Wtime2, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncz, Ztime2.st, ncw) / n
 
  para <- List2Vec(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  PLs <- matrix(0, len, len)
  for (i in 1:len) {
    for (j in i:len) {
      para1 <- para
      para1[i] <- para[i] + delta
      para1[j] <- para1[j] + delta
      result <- LambGeneric(para1, lamb.init, tol, iter, ncz, ncx, ncw, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2)
      para1.list <- Vec2List(para1, ncx, ncz, ncw)
      theta.input1 <- list(beta = para1.list$beta, phi = para1.list$phi, alpha = para1.list$alpha, 
                           Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
      PLs[i, j] <- LHGeneric(theta.input1, n, Z.st, Y.st, X.st, b, Ztime,  nk, Index0, model, Wtime, Xtime, Wtime2, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncz, Ztime2.st, ncw) / n
    }
  }
  
  pls <- rep(0, len)
  for (i in 1:len) {
    para1 <- para
    para1[i] <- para[i] + delta
    result <- LambGeneric(para1, lamb.init, tol, iter, ncz, ncx, ncw, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2)
    para1.list <- Vec2List(para1, ncx, ncz, ncw)
    theta.input1 <- list(beta = para1.list$beta, phi = para1.list$phi, alpha = para1.list$alpha, 
                         Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
    pls[i] <- LHGeneric(theta.input1, n, Z.st, Y.st, X.st, b, Ztime, nk, Index0, model, Wtime, Xtime, Wtime2, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncz, Ztime2.st, ncw) / n
  }
  
  I <- matrix(0, len, len)
  for (i in 1:len) {
    for (j in i:len) {
      I[i, j] <- - (PLs[i, j] - pls[i] - pls[j] + pl) / (delta ^ 2)
    }
  }
  I <- I + t(I) - diag(diag(I)) 
  V <- solve(I) / n;
  Valpha.name <- if (model == 1) paste("alpha:", alpha.name, sep = "") else "alpha"
  Vnames <- c(paste(rep("beta:", ncx), beta.names, sep = ""), paste(rep("phi:", ncw), phi.names, sep = ""),
              Valpha.name, "sigma.e", paste("Bsigma.", 1:p, sep = ""))
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
