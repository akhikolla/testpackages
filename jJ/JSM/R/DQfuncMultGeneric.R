
#=============== The DQ Function for Model I of Multiplicative Joint Model ===============#

DQfuncMultGeneric <- function (ptheta, theta, B.st, n, Y.st, b, model, Btime, Btime2,  Index, Index0, Ztime, Ztime2, nknot, nk, Index1, rho, d, wGQ, ncz, ncb, N, Y, B, ID, Index2 ) { # ptheta means "theta prime"
  
  pgamma <- ptheta$gamma
  pphi <- ptheta$phi
  palpha <- ptheta$alpha
  pYsigma2  <- (ptheta$Ysigma) ^ 2
  pBsigma2 <- (ptheta$Bsigma) ^ 2
  plamb <- ptheta$lamb
  gamma <- theta$gamma
  phi <- theta$phi
  alpha <- theta$alpha
  Ysigma2  <- (theta$Ysigma) ^ 2
  Bsigma2 <- (theta$Bsigma) ^ 2
  lamb <- theta$lamb
  
  M <- length(Index)
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))
  VY <- lapply(1:n, function(i) calc_VY(BTg[[i]], Bsigma2, Ysigma2)) 
  VB <-  lapply(1:n, function(i) calc_VB(M1 = Bsigma2,M2 =  BTg[[i]],  VY[[i]]))
  muB <-lapply(1:n, function(i) calc_muBMult(  Bsigma2,VY[[i]],BTg[[i]],Y.st[[i]] )+1 )
  bi.st <- lapply(1:n, function(i) calc_bi_st(v0=muB[[i]], b ,M = VB[[i]]) ) 

  bi <- do.call(rbind, bi.st) # n*nknot matrix #
  if( model==1) { 
    Btime.b <- as.vector(Btime %*% gamma) * bi # n*nknot matrix #
    Btime2.b <- as.vector(Btime2 %*% gamma) * bi[Index, ] # M*nknot matrix #
  }
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  
  Ztime_phi <- if (ncz > 0) Ztime %*% phi else rep(0, n)
  Ztime2_phi <- if (ncz > 0) Ztime2 %*% phi else rep(0, M)
  Ztime2_pphi <- if (ncz > 0) Ztime2 %*% pphi else rep(0, M)

  if( model==1) { 
    log.density1 <- log.lamb + as.vector(Ztime_phi) + alpha * Btime.b # n*nknot matrix #
    eta.s <- as.vector(Ztime2_phi) + alpha * Btime2.b # M*nknot matrix #
  } else if(model ==2){
   log.density1 <- log.lamb + as.vector(Ztime_phi) + alpha * bi # n*nknot matrix #
   eta.s <- as.vector(Ztime2_phi) + alpha * bi[Index, ] # M*nknot matrix #
  } else {
    stop("Invalid model type")
  }

  calc_expM2(eta.s)
  const <- matrix(0, n, nknot) # n*nknot matrix #
  const[nk != 0, ] <- calc_rowsum_mult((Index), lamb[Index1], eta.s )
  log.density2 <- - log(1 + rho * const) # n*nknot matrix # 
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*nknot matrix #
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*nknot matrix #
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*nknot matrix #
  
  len <- ncz + ncb + 3
  Q <- rep(0, len)
  
  Q[len] <- - n / sqrt(pBsigma2) + sum((Integral * (bi - 1) ^ 2) %*% wGQ) / (pBsigma2 ^ (3 / 2))
  Q[len - 1] <- - N / sqrt(pYsigma2) + sum(((Y - as.vector(B %*% pgamma) * bi[ID, ]) ^ 2 * Integral[ID, ]) %*% wGQ) / (pYsigma2 ^ (3 / 2))
  
  if( model==1) { 
    pBtime2.b <- as.vector(Btime2 %*% pgamma) * bi[Index, ] # M*nknot matrix #
    temp0 <- as.vector(Ztime2_pphi) + palpha * pBtime2.b # M*nknot matrix #
  } else { 
    temp0 <- as.vector(Ztime2_pphi) + palpha * bi[Index, ] # M*nknot matrix #
  }
  calc_expM2(temp0)

  calc_M1_M2_M3_Hadamard(temp0, CondExp, Integral,as.integer(Index-1))
  temp1 <- as.vector( temp0 %*% wGQ) # vector of length M #
  if( model==1) { 
    temp2 <- as.vector((pBtime2.b * temp0) %*% wGQ) # vector of length M # 
    calc_M1_a_M2_Hadamard( temp0, bi, palpha, as.integer(Index-1))
    temp3 <- calc_M1_M2_Hadamard_a(temp0,Btime2,wGQ, a = ncb)
  } else { 
    temp2 <- as.vector((bi[Index, ] * temp0) %*% wGQ) # vector of length M #
  }
  
  post1 <- calc_tapply_vect_sum(temp1,  as.integer(Index1-1))
  post2 <- calc_tapply_vect_sum(temp2,  as.integer(Index1-1))
  
  if (ncz > 0) {
    temp4 <- Ztime2 * temp1 # M*ncz matrix #
    post4 <- as.matrix(apply(temp4, 2, function(x) calc_tapply_vect_sum(x,  as.integer(Index1-1)))) # n_u*ncz matrix #
    Q[(ncb + 1) : (ncb + ncz)] <- colSums(d * Ztime) - colSums(Index2 * post4 / post1) # vector of length ncz #
  }
  
  temp5 <- lapply(1:ncb, function(i) (Y - as.vector(B %*% pgamma) * bi[ID, ]) * B[, i] * bi[ID, ]) # N*nknot matrices #
  
  if( model==1) { 
    post3 <- as.matrix(apply(temp3, 2, function(x) calc_tapply_vect_sum(x,  as.integer(Index1-1)))) # n_u*ncb matrix # 
  }

  post5 <- unlist(lapply(temp5, function(x) sum((x * Integral[ID, ]) %*% wGQ))) # vector of length ncb #
  post.bi <- as.vector((Integral * bi) %*% wGQ) # vector of length n #
  
  if( model==1) { 
    Q[1 : ncb] <- palpha * colSums(d * post.bi * Btime) - colSums(Index2 * post3 / post1) + post5 / pYsigma2
    Q[ncb + ncz + 1] <- sum(d * post.bi * as.vector(Btime %*% pgamma)) - sum(Index2 * post2 / post1)
  } else {
    Q[1 : ncb] <- post5 / pYsigma2
    Q[ncb + ncz + 1] <- sum(d * post.bi) - sum(Index2 * post2 / post1)
  }

  return(Q)
}
