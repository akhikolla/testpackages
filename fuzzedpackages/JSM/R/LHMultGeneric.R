
#=============== Function to Calculate the Likelihood Value for Model I ===============#
#================ Multiplicative Joint Modeling ===============#

LHMultGeneric <- function (theta, B.st, n, Y.st, b, model, Btime, Btime2, Index, Index0, Ztime, Ztime2, nknot, nk, Index1, rho, d, wGQ, ncz) {
  
  gamma <- theta$gamma
  phi <- theta$phi
  alpha <- theta$alpha
  Ysigma2 <- (theta$Ysigma) ^ 2
  Bsigma2 <- (theta$Bsigma) ^ 2
  lamb <- theta$lamb
  
  M <- length(Index)
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))
  VY <- lapply(1:n, function(i) calc_VY(BTg[[i]], Bsigma2, Ysigma2)) 
  VB <-  lapply(1:n, function(i) calc_VB(M1 = Bsigma2,M2 =  BTg[[i]],  VY[[i]]))
  muB <-lapply(1:n, function(i) calc_muBMult(  Bsigma2,VY[[i]],BTg[[i]],Y.st[[i]] )+1 )
  bi.st <- lapply(1:n, function(i) calc_bi_st(v0=muB[[i]], b ,M = VB[[i]]) ) 

  bi <- do.call(rbind, bi.st) # n*nknot matrix #
  if (model==1){
    Btime.b <- as.vector(Btime %*% gamma) * bi # n*nknot matrix #
    Btime2.b <- as.vector(Btime2 %*% gamma) * bi[Index, ] # M*nknot matrix #
  }
  
  Ztime_phi <- if (ncz > 0) Ztime %*% phi else rep(0, n)
  Ztime2_phi <- if (ncz > 0) Ztime2 %*% phi else rep(0, M)
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  if (model==1){
    log.density1 <- log.lamb + as.vector(Ztime_phi) + alpha * Btime.b # n*nknot matrix #
    exp.es <- alpha * Btime2.b + as.vector(Ztime2_phi)
  } else {
    log.density1 <- log.lamb + as.vector(Ztime_phi) + alpha * bi # n*nknot matrix # 
    exp.es <- alpha * bi[Index, ] + as.vector(Ztime2_phi)
  }
  calc_expM2(exp.es);

  const <- matrix(0, n, nknot) # n*nknot matrix #
  const[nk != 0, ] <- calc_rowsum( (Index),  exp.es * lamb[Index1])

  log.density2 <- - log(1 + rho * const) # n*nknot matrix # 
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*nknot matrix #
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  
  f.long <- sapply(1:n, function(i) calc_MVND(Y.st[[i]], as.vector(BTg[[i]]), VY[[i]]))
  return(sum(log(f.long * deno / sqrt(pi))))
}
