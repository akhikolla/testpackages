
#========== Function to Obtain Lamb Given Other Finite Dimensional Parameters ==========#
#=============== Model I for Multiplicative Joint Modeling ===============#

LambMultGeneric <- function (para, lamb.init, tol, iter, ncz, ncb, B.st, n, Y.st, b, model, Btime, Btime2, Index, Ztime, Ztime2, Index0, nknot, nk, Index1, rho, d, wGQ, Index2){
  
  para.list <- Vec2ListMult(para, ncz, ncb)
  gamma <- para.list$gamma
  phi <- para.list$phi
  alpha <- para.list$alpha
  Ysigma2 <- (para.list$Ysigma) ^ 2
  Bsigma2 <- (para.list$Bsigma) ^ 2
  
  M <- length(Index)
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))
  VY <- lapply(1:n, function(i) calc_VY(BTg[[i]], Bsigma2, Ysigma2)) 
  VB <-  lapply(1:n, function(i) calc_VB(M1 = Bsigma2,M2 =  BTg[[i]],  VY[[i]]))
  muB <-lapply(1:n, function(i) calc_muBMult(  Bsigma2,VY[[i]],BTg[[i]],Y.st[[i]] )+1 )
  bi.st <- lapply(1:n, function(i) calc_bi_st(v0=muB[[i]], b ,M = VB[[i]]) ) 

  bi <- do.call(rbind, bi.st) # n*nknot matrix #

  Ztime_phi <- if (ncz > 0) Ztime %*% phi else rep(0, n)
  Ztime2_phi <- if (ncz > 0) Ztime2 %*% phi else rep(0, M)
  
  if (model == 1){
    Btime.b <- as.vector(Btime %*% gamma) * bi # n*nknot matrix #
    Btime2.b <- as.vector(Btime2 %*% gamma) * bi[Index, ] # M*nknot matrix #
  
    eta.h <- as.vector(Ztime_phi) + alpha * Btime.b # n*nknot matrix #
    exp.es <- exp( as.vector(Ztime2_phi) + alpha * Btime2.b) # M*nknot matrix #
  } else if(model ==2){ 
    eta.h <- as.vector(Ztime_phi) + alpha * bi # n*nknot matrix #
    exp.es <- exp(as.vector(Ztime2_phi) + alpha * bi[Index, ]) # M*nknot matrix #
  } else {
    stop("Invalid model type")
  }   
  
  lamb.old <- lamb.init
  err <- 1
  
  for (step in 1:iter) {
    Lamb.old <- cumsum(lamb.old)
    
    log.lamb <- log(lamb.old[Index0])
    log.lamb[is.na(log.lamb)] <- 0
    log.density1 <- log.lamb + eta.h # n*nknot matrix #
    const <- matrix(0, n, nknot) # n*nknot matrix # 
    const[nk != 0, ] <- calc_rowsum_mult((Index), lamb.old[Index1], exp.es)
    log.density2 <- - log(1 + rho * const) # n*nknot matrix # 
    log.survival <- if(rho > 0) - log(1 + rho * const) / rho else - const # n*nknot matrix # 
    
    f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
    deno <- as.vector(f.surv %*% wGQ) # vector of length n #
    Integral <- f.surv / deno # n*nknot matrix #
    CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
     
    tempLamb0 <- exp.es; tempLamb0[1] = tempLamb0[1] + 0 # "touch the variable"
    calc_M1_M2_M3_Hadamard(tempLamb0, CondExp , Integral, as.integer(Index - 1))
    tempLamb <- calc_M_v(v = wGQ, M = tempLamb0)   
    postLamb <- calc_tapply_vect_sum(v1 = tempLamb, v2 = as.integer(Index1 - 1))
    lamb.new <- Index2 / postLamb
    
    Lamb.new <- cumsum(lamb.new)
    err <- max(abs((Lamb.new - Lamb.old) / Lamb.old))
    if (step > 3 & err < tol) break
    lamb.old <- lamb.new
  }
  converge <- as.numeric(step <= iter & err < tol)
  return(list(lamb = lamb.new, converge = converge))
}
