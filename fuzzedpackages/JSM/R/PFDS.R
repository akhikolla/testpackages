
#========== Differentiate the S function with Forward Difference for Model I & II ==========#

PFDS <- function (model, theta, cvals, ncx, ncz, ncw, p, varNames, n, Z.st, Y.st, X.st, Ztime, nk, Wtime, Wtime2, Xtime, Xtime2, GQ,  rho, d, wGQ, ncz2, b, Ztime2.st,  X, Y, ID, N, Indcs, Z) {

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

  S <- Sfunc(model, theta, n, Z.st, Y.st, X.st, Ztime, nk, Wtime,  Wtime2, Xtime, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncx, ncw, p, ncz, ncz2, b, Ztime2.st, Index0, X, Y, ID, N, Index2, Z)
  
  para <- List2Vec(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  #===== Calculate the derivative of the score vector using forward difference =====#
  DS <- matrix(0, len, len)
  
  for (i in 1:len) {
    para1 <- para
    para1[i] <- para[i] + delta
    result <- LambGeneric(para1, lamb.init, tol, iter, ncz, ncx, ncw, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2)
    para1.list <- Vec2List(para1, ncx, ncz, ncw)
    theta.input1 <- list(beta = para1.list$beta, phi = para1.list$phi, alpha = para1.list$alpha, 
                         Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
    S1 <- Sfunc(model, theta.input1, n, Z.st, Y.st, X.st, Ztime, nk, Wtime, Wtime2, Xtime, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncx, ncw, p, ncz, ncz2, b, Ztime2.st, Index0, X, Y, ID, N, Index2, Z)
    DS[i, ] <- (S1 - S) / delta
  }
  
  #========== make the DS matrix symmetric ==========#
  DS <- (DS + t(DS)) / 2 
  V <- -solve(DS);
  Valpha.name <- if (model == 1) paste("alpha:", alpha.name, sep = "") else "alpha"
  Vnames <- c(paste(rep("beta:", ncx), beta.names, sep = ""), paste(rep("phi:", ncw), phi.names, sep = ""),
              Valpha.name, "sigma.e", paste("Bsigma.", 1:p, sep = ""))
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
