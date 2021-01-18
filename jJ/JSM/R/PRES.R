
#========== Differentiate the S function with Richardson Extrapolation for Model I & II ==========#

PRES <- function (model, theta, cvals, ncz, ncx, ncw,n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, wGQ, GQ, d, p, ncz2, X, Y, ID, N, Z, varNames, Indcs) {     
 
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

  para <- List2Vec(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  #===== Calculate the derivative of the score vector using forward difference =====#
  DS <- matrix(0, len, len)
  
  for (i in 1:len) {
    para1 <- para2 <- para3 <- para4 <- para
    para1[i] <- para[i] - 2 * delta
    para2[i] <- para[i] - delta
    para3[i] <- para[i] + delta
    para4[i] <- para[i] + 2 * delta
    result1 <- LambGeneric(para1, lamb.init, tol, iter, ncz, ncx, ncw, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2) 
    result2 <- LambGeneric(para2, lamb.init, tol, iter, ncz, ncx, ncw, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2)
    result3 <- LambGeneric(para3, lamb.init, tol, iter, ncz, ncx, ncw, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2)
    result4 <- LambGeneric(para4, lamb.init, tol, iter, ncz, ncx, ncw, n, Z.st, Y.st, X.st, b, Ztime, Ztime2.st, nk, Wtime, Xtime, Wtime2, Xtime2, rho, Index0, Index1, Index, wGQ, model, GQ, d, Index2)
    list1 <- Vec2List(para1, ncx, ncz, ncw)
    list2 <- Vec2List(para2, ncx, ncz, ncw)
    list3 <- Vec2List(para3, ncx, ncz, ncw)
    list4 <- Vec2List(para4, ncx, ncz, ncw)
    theta.input1 <- list(beta = list1$beta, phi = list1$phi, alpha = list1$alpha, Ysigma = list1$Ysigma, Bsigma = list1$Bsigma, lamb = result1$lamb)
    theta.input2 <- list(beta = list2$beta, phi = list2$phi, alpha = list2$alpha, Ysigma = list2$Ysigma, Bsigma = list2$Bsigma, lamb = result2$lamb)
    theta.input3 <- list(beta = list3$beta, phi = list3$phi, alpha = list3$alpha, Ysigma = list3$Ysigma, Bsigma = list3$Bsigma, lamb = result3$lamb)
    theta.input4 <- list(beta = list4$beta, phi = list4$phi, alpha = list4$alpha, Ysigma = list4$Ysigma, Bsigma = list4$Bsigma, lamb = result4$lamb)
    S1 <- DQfuncGeneric(model = model, ptheta = theta.input1, theta = theta.input1, n =  n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, Ztime = Ztime, nk = nk, Wtime = Wtime, Wtime2 = Wtime2, Xtime = Xtime, Xtime2 = Xtime2, GQ = GQ, Index = Index, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncx = ncx, ncw = ncw, p = p, ncz = ncz, ncz2 = ncz2, b = b, Ztime2.st = Ztime2.st, Index0 = Index0, X = X, Y = Y , ID = ID, N = N, Index2 = Index2, Z= Z)
    S2 <- DQfuncGeneric(model = model, ptheta = theta.input2, theta = theta.input2, n =  n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, Ztime = Ztime, nk = nk, Wtime = Wtime, Wtime2 = Wtime2, Xtime = Xtime, Xtime2 = Xtime2, GQ = GQ, Index = Index, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncx = ncx, ncw = ncw, p = p, ncz = ncz, ncz2 = ncz2, b = b, Ztime2.st = Ztime2.st, Index0 = Index0, X = X, Y = Y , ID = ID, N = N, Index2 = Index2, Z= Z)
    S3 <- DQfuncGeneric(model = model, ptheta = theta.input3, theta = theta.input3, n =  n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, Ztime = Ztime, nk = nk, Wtime = Wtime, Wtime2 = Wtime2, Xtime = Xtime, Xtime2 = Xtime2, GQ = GQ, Index = Index, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncx = ncx, ncw = ncw, p = p, ncz = ncz, ncz2 = ncz2, b = b, Ztime2.st = Ztime2.st, Index0 = Index0, X = X, Y = Y , ID = ID, N = N, Index2 = Index2, Z= Z)
    S4 <- DQfuncGeneric(model = model, ptheta = theta.input4, theta = theta.input4, n =  n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, Ztime = Ztime, nk = nk, Wtime = Wtime, Wtime2 = Wtime2, Xtime = Xtime, Xtime2 = Xtime2, GQ = GQ, Index = Index, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncx = ncx, ncw = ncw, p = p, ncz = ncz, ncz2 = ncz2, b = b, Ztime2.st = Ztime2.st, Index0 = Index0, X = X, Y = Y , ID = ID, N = N, Index2 = Index2, Z= Z)
    DS[i, ] <- (S1 - 8 * S2 + 8 * S3 - S4) / (12 * delta)
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
