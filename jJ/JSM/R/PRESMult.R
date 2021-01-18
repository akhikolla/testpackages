
#========== Differentiate the S function with Richardson Extrapolation for Model I & II ==========#
#========== Multipicative Joint Modeling ==========#

PRESMult <- function (model, theta, tol, iter, delta, ncz, ncb, B.st, n, Y.st, b, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2, alpha.name, phi.names, N, Y, B, ID) {     

para <- List2VecMult(theta)
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
  result1 <- LambMultGeneric(para1, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st, b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2 )
  result2 <- LambMultGeneric(para2, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st, b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2 )
  result3 <- LambMultGeneric(para3, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st, b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2 )
  result4 <- LambMultGeneric(para4, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st, b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2 )
  list1 <- Vec2ListMult(para1, ncz, ncb)
  list2 <- Vec2ListMult(para2, ncz, ncb)
  list3 <- Vec2ListMult(para3, ncz, ncb)
  list4 <- Vec2ListMult(para4, ncz, ncb)
  theta.input1 <- list(gamma = list1$gamma, phi = list1$phi, alpha = list1$alpha, 
                       Ysigma = list1$Ysigma, Bsigma = list1$Bsigma, lamb = result1$lamb)
  theta.input2 <- list(gamma = list2$gamma, phi = list2$phi, alpha = list2$alpha, 
                       Ysigma = list2$Ysigma, Bsigma = list2$Bsigma, lamb = result2$lamb)
  theta.input3 <- list(gamma = list3$gamma, phi = list3$phi, alpha = list3$alpha, 
                       Ysigma = list3$Ysigma, Bsigma = list3$Bsigma, lamb = result3$lamb)
  theta.input4 <- list(gamma = list4$gamma, phi = list4$phi, alpha = list4$alpha, 
                       Ysigma = list4$Ysigma, Bsigma = list4$Bsigma, lamb = result4$lamb)
  S1 <- SfuncMult(model, theta.input1, B.st = B.st, n = n, Y.st = Y.st, b = b, Btime = Btime, Btime2 =  Btime2, Index =  Index, Index0 =  Index0, Ztime =  Ztime, Ztime2 = Ztime2, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncz = ncz, ncb = ncb, N = N, Y = Y, B = B, ID = ID, Index2 = Index2)
  S2 <- SfuncMult(model, theta.input2, B.st = B.st, n = n, Y.st = Y.st, b = b, Btime = Btime, Btime2 =  Btime2, Index =  Index, Index0 =  Index0, Ztime =  Ztime, Ztime2 = Ztime2, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncz = ncz, ncb = ncb, N = N, Y = Y, B = B, ID = ID, Index2 = Index2)
  S3 <- SfuncMult(model, theta.input3, B.st = B.st, n = n, Y.st = Y.st, b = b, Btime = Btime, Btime2 =  Btime2, Index =  Index, Index0 =  Index0, Ztime =  Ztime, Ztime2 = Ztime2, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncz = ncz, ncb = ncb, N = N, Y = Y, B = B, ID = ID, Index2 = Index2)
  S4 <- SfuncMult(model, theta.input4, B.st = B.st, n = n, Y.st = Y.st, b = b, Btime = Btime, Btime2 =  Btime2, Index =  Index, Index0 =  Index0, Ztime =  Ztime, Ztime2 = Ztime2, nknot = nknot, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncz = ncz, ncb = ncb, N = N, Y = Y, B = B, ID = ID, Index2 = Index2)
  DS[i, ] <- (S1 - 8 * S2 + 8 * S3 - S4) / (12 * delta)
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
