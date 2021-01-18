
#=============== EM iteration Using Adaptive Gaussian Quadrature for Model I with NMRE ===============#
#=============== Transformation model is fitted for the survival part ===============#

EMiterMultGeneric <- function (theta.old, B.st, n, Y.st, b, model, Btime, Btime2, Index, Index0, Ztime, Ztime2, nknot, nk, Index1, rho, d, wGQ, ID, ncb, B, Y, N, ncz, Ztime22, Index2, B2, Btime22) { # Use apply instead of matrix calculation #
  
  # Get Old Estimates #
  gamma <- theta.old$gamma
  phi.old <- theta.old$phi
  alpha.old <- theta.old$alpha
  Ysigma2.old <- (theta.old$Ysigma) ^ 2
  Bsigma2.old <- (theta.old$Bsigma) ^ 2
  lamb.old <- theta.old$lamb
  
  M <- length(Index)
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma)) 
  VY <- lapply(1 : n, function(i) calc_VY( BTg[[i]], Bsigma2.old, Ysigma2.old) )
  VB <-  lapply(1 : n, function(i) calc_VB(M1 = Bsigma2.old, M2 = BTg[[i]], M3 = VY[[i]]))
  muB <- lapply(1 : n, function(i) calc_muBMult(Bsigma2.old, VY[[i]], BTg[[i]], Y.st[[i]]) + 1)
  bi.st <- lapply(1 : n, function(i) calc_bi_st(v0 = muB[[i]], b, M = VB[[i]]) ) 
 
  bi <- do.call(rbind, bi.st) # n*nknot matrix #
  if (model == 1) {
    Btime.b <- as.vector(Btime %*% gamma) * bi # n*nknot matrix #
    Btime2.b <- as.vector(Btime2 %*% gamma) * bi[Index, ] # M*nknot matrix #
  } else if (model == 2){
    # Do nothing
  } else {
    stop("Invalid model type")
  }    
  log.lamb <- log(lamb.old[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  
  Ztime_phi.old <- if (ncz > 0) Ztime %*% phi.old else rep(0, n)
  Ztime2_phi.old <- if (ncz > 0) Ztime2 %*% phi.old else rep(0, M)
  
  if (model == 1){ 
    log.density1 <- log.lamb + as.vector(Ztime_phi.old) + alpha.old * Btime.b # n*nknot matrix #
    exp.es <- as.vector(Ztime2_phi.old) + alpha.old * Btime2.b # M*nknot matrix #
  } else { 
    log.density1 <- log.lamb + as.vector(Ztime_phi.old) + alpha.old * bi # n*nknot matrix #
    exp.es <- as.vector(Ztime2_phi.old) + alpha.old * bi[Index, ] # M*nknot matrix #
  }
  # exp.es <- exp(eta.s) # M*nknot matrix #
  calc_expM2(exp.es)

  const <- matrix(0, n, nknot) # n*nknot matrix #
  const[nk != 0, ] <- calc_rowsum_mult((Index), lamb.old[Index1], exp.es) # n*nknot matrix #  
  log.density2 <- -log(1 + rho * const) # n*GQ matrix # 
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*nknot matrix #
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*nknot matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*nknot matrix #
  
  f.long <- sapply(1 : n, function(i) calc_MVND(Y.st[[i]], as.vector(BTg[[i]]), VY[[i]]))
  lgLik <- sum(log(f.long * deno / sqrt(pi)))
  
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*nknot matrix #
  
  #========== Update Bsigma ==========#
  Bsigma2.new <- mean((Integral * (bi - 1) ^ 2) %*% wGQ) 

  if (model == 2) {
    #========== Update gamma: the linear regresion coefficents of regression Yi on E(bi)*B ==========#
    post.bi <- as.vector((Integral * bi) %*% wGQ) # vector of length n #
    post.bi2 <- as.vector((Integral * bi ^ 2) %*% wGQ) # vector of length n #
    gamma <- as.vector(solve(matrix(colSums(post.bi2[ID] * B2), ncb)) %*% colSums(post.bi[ID] * B * Y)) 
    # vector of length ncb #
  }
  
  #========== Update Ysigma ==========#
  post.resid <- ((Y - as.vector(B %*% gamma) * bi[ID, ]) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Ysigma2.new <- sum(post.resid) / N
  
  #========== calculate the score and gradient of phi and alpha ==========# 
  CondExp2 <- CondExp[nk != 0, ]
  temp0 <- exp.es * lamb.old[Index1] 

  if (model == 2) { 
    temp0c <- bi[Index, ]* temp0
    temp2 <- CondExp2 * calc_rowsum(temp0c,  v = Index) # n*nknot matrix # 
    temp4 <- calc_mult_rowsum2(Index, bi[Index, ], temp0c, A = CondExp2)
  } else {
    temp0c <- Btime2.b * temp0;
    temp2 <- CondExp2 * calc_rowsum(temp0c,  v = Index) # n*nknot matrix # 
    temp4 <- calc_mult_rowsum2(A = CondExp2, v = Index, Btime2.b, temp0c) 
  }  

  Integral2 <- Integral[nk != 0, ]
  post2 <- sum((temp2 * Integral2) %*% wGQ)
  post4 <- sum((temp4 * Integral2) %*% wGQ)

  if (ncz > 0) {
    temp1 <- lapply(1 : ncz, function(i)  calc_mult_rowsum1((Index), Ztime2[, i] , temp0, CondExp2)) # n*nknot matrices # 
    temp3 <- lapply(1 : (ncz ^ 2), function(i) calc_mult_rowsum1(Index, Ztime22[, i], temp0, A = CondExp2)) # n*nknot matrices # 
    temp5 <- lapply(1 : (ncz), function(i) calc_mult_rowsum1(Index, Ztime2[, i], temp0c, A = CondExp2)) # n*nknot matrices #
    post1 <- unlist(lapply(temp1, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncz #
    post3 <- unlist(lapply(temp3, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncz^2 #
    post5 <- unlist(lapply(temp5, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncz #
  }
  
  if (model == 1) {
    post.bi <- as.vector((Integral * bi) %*% wGQ) # vector of length n #
  }
  
  if (model == 1) {
    alphaScore <- sum(d * post.bi * as.vector(Btime %*% gamma)) - post2
  } else {
    alphaScore <- sum(d * post.bi) - post2
  }
  
  if (ncz > 0) {
    phiScore <- colSums(d * Ztime) - post1 # vector of length ncz #
    pa.score <- c(phiScore, alphaScore)
    pa.info <- matrix(0, (ncz + 1), (ncz + 1)) # (ncz+1)*(ncz+1) matrix #
    pa.info[1 : ncz, 1 : ncz] <- - post3
    pa.info[(ncz + 1), (ncz + 1)] <- - post4
    pa.info[(ncz + 1), 1:ncz] <- - post5
    pa.info[1 : ncz, (ncz + 1)] <- - post5
    
    #=============== Update phi and alpha ===============#
    pa.old <- c(phi.old, alpha.old) # vector of length (ncz+1) #
    paSVD <- svd(pa.info)
    pa.info.inv <- paSVD$v %*% diag(1 / paSVD$d) %*% t(paSVD$u)
    pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncz+1) #
    phi.new <- pa.new[1 : ncz] 
    alpha.new <- pa.new[ncz + 1]
  } else {
    alpha.new <- alpha.old - alphaScore / (-post4)
    phi.new <- phi.old
  }

  Ztime2_phi.new <- if (ncz > 0) Ztime2 %*% phi.new else rep(0, M)
  if ( model == 1) {
    #========== calculate the score of gamma ==========#
    exp.es.n1 <- as.vector(Ztime2_phi.new) + alpha.new * Btime2.b # M*nknot matrix #
    #exp.es.n1 <- exp(eta.s.n1) # M*nknot matrix #
    calc_expM2(exp.es.n1)
    temp0b <- exp.es.n1 * lamb.old[Index1] * alpha.new * bi[Index, ]
    temp6 <- calc_mult_rowsum3( Index, Btime2, temp0b, A = CondExp2, ncb)  
    temp0c <-  alpha.new * bi[Index, ] *temp0b
    temp7 <- calc_mult_rowsum3( Index, Btime22, temp0c, A = CondExp2, ncb ^ 2)
    temp8 <- lapply(1:ncb, function(i) (Y - as.vector(B %*% gamma) * bi[ID, ]) * B[, i] * bi[ID, ]) 
    # N*nknot matrices #
    post6 <- unlist(lapply(temp6, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncb #
    post7 <- unlist(lapply(temp7, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncb^2 #
    post8 <- unlist(lapply(temp8, function(x) sum((x * Integral[ID,]) %*% wGQ))) # vector of length ncb #
    post.bi2 <- as.vector((Integral * bi ^ 2) %*% wGQ) # vector of length n #
  
    gammaScore <- post8 / Ysigma2.new + alpha.new * colSums(d * post.bi * Btime) - post6 # vector of length ncb #
    gammaInfo <- - colSums(post.bi2[ID] * B2) / Ysigma2.new - post7
    gammaInfo <- matrix(gammaInfo, nrow = ncb)
  
    #========== Update gamma ==========#
    gammaSVD <- svd(gammaInfo)
    gammaInfo.inv <- gammaSVD$v %*% diag(1 / gammaSVD$d) %*% t(gammaSVD$u)
    gamma <- as.vector(gamma - gammaInfo.inv %*% gammaScore) # vector of length ncb #
  } 

  #========== Calculate the new lambda with new parameters ==========# 
  if (model == 1) {
    Btime2.bnew <- as.vector(Btime2 %*% gamma) * bi[Index, ] # M*nknot matrix #
    eta.s.n2 <- as.vector(Ztime2_phi.new) + alpha.new * Btime2.bnew # M*nknot matrix #
  } else {
    eta.s.n2 <- as.vector(Ztime2_phi.new) + alpha.new * bi[Index, ] # M*nknot matrix #
  }

  calc_expM2(eta.s.n2)
  calc_M1_M2_M3_Hadamard(eta.s.n2, CondExp ,  Integral, as.integer(Index - 1))
  tempLamb <- calc_M_v(v = wGQ, M = eta.s.n2) 
  postLamb <- calc_tapply_vect_sum(v1 = tempLamb, v2 = as.integer(Index1 - 1)) # vector of length n_u #
  lamb.new <- Index2 / postLamb
  
  result <- list(gamma = gamma, phi = phi.new, alpha = alpha.new, Ysigma = sqrt(Ysigma2.new), 
                 Bsigma = sqrt(Bsigma2.new), lamb = lamb.new, lgLik = lgLik, est.bi = post.bi)
  return(result)
}
