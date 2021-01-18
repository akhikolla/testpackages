
#========== Differentiate the S function with Forward Difference for Model I & II ==========#
#========== Multipicative Joint Modeling ==========#

PFDSMult <- function (model, theta, tol, iter, delta, ncz, ncb, B.st, n, Y.st, b, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2, alpha.name, phi.names, N, Y, B, ID) {
  
  S <- SfuncMult(model, theta, B.st = B.st, n = n, Y.st = Y.st, b = b, Btime = Btime, Btime2 =  Btime2, Index =  Index, Index0 =  Index0, Ztime =  Ztime, Ztime2 = Ztime2, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncz = ncz, ncb = ncb, N = N, Y = Y, B = B, ID = ID, Index2 = Index2)
  
  para <- List2VecMult(theta)
  lamb.init <- theta$lamb
  len <- length(para)
  
  #===== Calculate the derivative of the score vector using forward difference =====#
  DS <- matrix(0, len, len)
  
  for (i in 1:len) {
    para1 <- para
    para1[i] <- para[i] + delta
    result <-   LambMultGeneric(para1, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st, b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2)
    para1.list <- Vec2ListMult(para1, ncz, ncb)
    theta.input1 <- list(gamma = para1.list$gamma, phi = para1.list$phi, alpha = para1.list$alpha, 
                         Ysigma = para1.list$Ysigma, Bsigma = para1.list$Bsigma, lamb = result$lamb)
    S1 <- SfuncMult(model, theta.input1, B.st = B.st, n = n, Y.st = Y.st, b = b, Btime = Btime, Btime2 =  Btime2, Index =  Index, Index0 =  Index0, Ztime =  Ztime, Ztime2 = Ztime2, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncz = ncz, ncb = ncb, N = N, Y = Y, B = B, ID = ID, Index2 = Index2)
    DS[i, ] <- (S1 - S) / delta
  }
  
  #========== make the DS matrix symmetric ==========#
  DS <- (DS + t(DS)) / 2
  V <- -solve(DS);
  Valpha.name <- if (model == 1) paste("alpha:", alpha.name, sep = "") else "alpha"
  Vnames <- c(paste("gamma.", 1:ncb, sep = ""), paste(rep("phi:", ncz), phi.names, sep = ""), Valpha.name, 
              "sigma.e", "sigma.b")
  dimnames(V) <- list(Vnames, Vnames)
  
  return(V)
}
