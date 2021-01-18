
#=============== EM iteration Using Adaptive Gaussian Quadrature for Model I&II ===============#
#=============== Transformation model is fitted for the survival part ===============#

EMiterTMGeneric <- function (theta.old, n, Z.st, Ztime, Ztime2.st, nk, Wtime2, Xtime2, GQ, rho, wGQ, d, Y.st, X.st, ncz, ncz2, b, model, Wtime, Xtime, X, Y, ID, N, ncw, Wtime22, ncx, Xtime22, Z, X2.sum, Indcs){ # Use apply instead of matrix calculation #

    # Get Old Estimates #
  beta.old <- theta.old$beta
  Ysigma2.old <- (theta.old$Ysigma) ^ 2
  Bsigma.old <- theta.old$Bsigma
  phi.old <- theta.old$phi
  alpha.old <- theta.old$alpha

  Index = Indcs$Index  
  Index0 = Indcs$Index0
  Index1 = Indcs$Index1
  Index2 = Indcs$Index2

  lamb.old <- theta.old$lamb
  
  M <- nrow(Xtime2)
  
  VY <- lapply(1 : n, function(i) calc_VY(M = Z.st[[i]], A = Bsigma.old, b = Ysigma2.old))  
  VB <-  lapply(1 : n, function(i) calc_VB(Bsigma.old, M2 =  Z.st[[i]], M3 = VY[[i]])) 
  muB <- lapply(1 : n, function(i) calc_muB(BSold = Bsigma.old, Zst = Z.st[[i]], Yst = Y.st[[i]], betaold = beta.old, VY = VY[[i]], Xst = X.st[[i]]))
  bi.st <- lapply(1 : n, function(i) calc_bi_st(v0 = muB[[i]], b, M = VB[[i]]) ) 
 
  bi <- do.call(rbind, bi.st)
  Ztime.b <- do.call(rbind, lapply(1 : n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix #
  Ztime2.b <-fast_lapply_length(Ztime2.st, bi.st, (1 : n)[nk != 0] - 1)# M*GQ matrix #
  
  log.lamb <- log(lamb.old[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  
  Wtime_phi.old <- if (ncw > 0) Wtime %*% phi.old else rep(0, n)
  Wtime2_phi.old <- if (ncw > 0) Wtime2 %*% phi.old else rep(0, M)

  if (model == 2) {
    log.density1 <- log.lamb + as.vector(Wtime_phi.old) + alpha.old * Ztime.b # n*GQ matrix #
    eta.s <- as.vector(Wtime2_phi.old) + alpha.old * Ztime2.b # M*GQ matrix #  
  } else if (model == 1) {
    log.density1 <- log.lamb + as.vector(Wtime_phi.old + alpha.old * Xtime %*% beta.old) + alpha.old * Ztime.b # n*GQ matrix #
    eta.s <- as.vector(Wtime2_phi.old + alpha.old * Xtime2 %*% beta.old) + alpha.old * Ztime2.b
  } else {
    stop("Invalid model type")
  } 

  const <- matrix(0, n, GQ) # n*GQ matrix #

  calc_expM2(eta.s)
  temp0a <- eta.s * lamb.old[Index1]

  const[nk != 0, ] <- calc_rowsum( (Index), temp0a)
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  log.survival <- if (rho > 0) log.density2 / rho else - const # n*GQ matrix # 
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*GQ matrix f(bi|Oi) #
  
  f.long <- sapply(1 : n, function(i) calc_MVND(Y.st[[i]], as.vector(X.st[[i]] %*% beta.old), VY[[i]]))

  lgLik <- sum(log(f.long * deno / (pi ^ (ncz / 2))))
  
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
  
  # post.bi <- Integral %*% (t(bi) * wGQ) # n*(n*ncz) matrix #
  post.bi <- calc_M1timesM2v(Integral, bi, wGQ)
  post.bi <- if(ncz > 1) {
		t(sapply(1 : n, function(i) post.bi[i, ((i - 1) * ncz + 1) : (i * ncz)])) 
		} else { 
		  matrix(diag(post.bi), nrow = n) # n*ncz matrix Ehat(bi) #
	}  
  #========== Update Bsigma ==========#
  if (ncz > 1) {
    tempB <-  fast_rbind_lapply_outerprod( bi.st )    # (n*ncz^2)*GQ matrix #      
  } else {
    tempB <- bi ^ 2
  }
  # post.bi2 <- Integral %*% (t(tempB) * wGQ) # n*(n*ncz^2) matrix #
  post.bi2 <- calc_M1timesM2v(Integral, tempB, wGQ)
  post.bi2 <- if(ncz > 1) {
		t(sapply(1 : n, function(i) post.bi2[i, ((i - 1) * ncz2 + 1) : (i * ncz2)])) 
	} else { 
	  matrix(diag(post.bi2), nrow = n) # n*(ncz^2) matrix Ehat(bibi^T) #
  }
  Bsigma.new <- if (ncz > 1) matrix(colMeans(post.bi2), ncz, ncz) else mean(post.bi2) # ncz*ncz matrix #
  
  if(model == 2) {
    #========== Update beta: the linear regresion coefficents of regression Yi-E(Zi*bi) on X_i ==========#
    tempX <- Y - rowSums(Z * post.bi[ID, ]) # vector of length N #
    beta.old <- as.vector(qr.solve(X, tempX)); # vector of length ncx #
    beta.new <- beta.old
  }
  
  #========== Update Ysigma ==========# 
  Ymu <- as.vector(X %*% beta.old) + do.call(rbind, lapply(1:n, function(i) Z.st[[i]] %*% bi.st[[i]])) # N*GQ matrix #

  post.resid <- ((Y - Ymu) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Ysigma2.new <- sum(post.resid) / N
  
  #========== calculate the score and gradient of phi and alpha ==========#
  CondExp2 <- CondExp[nk != 0, ]
  if (model == 2) { 
    temp0b <- Ztime2.b * temp0a 
    temp1 <- CondExp2 * calc_rowsum( (Index), temp0b) 
    temp2 <- calc_mult_rowsum2(v = Index, A = CondExp2, M = temp0b, L = Ztime2.b)
    if (ncw > 0) {
      temp3 <- lapply(1:(ncw), function(i) calc_mult_rowsum1(v = Index, u = Wtime2[, i], A = CondExp2,  M =temp0a))
      temp4 <- lapply(1:(ncw^2), function(i) calc_mult_rowsum1(v = Index, u = Wtime22[, i], A = CondExp2, M = temp0a)) 
    }
    temp0c <- Ztime2.b *temp0a
  } else {
    XZb2 <- as.vector(Xtime2 %*% beta.old) + Ztime2.b # M*GQ matrix # 
    temp0b <- XZb2 * temp0a; 
    temp1 <- CondExp2 * calc_rowsum( (Index), temp0b)
    temp2 <- calc_mult_rowsum2(v = Index, A = CondExp2, M = temp0b, XZb2)
    if (ncw > 0) {
      temp3 <- lapply(1:(ncw), function(i) calc_mult_rowsum1(v = Index, u = Wtime2[, i], A = CondExp2,  M = temp0a))
      temp4 <- lapply(1:(ncw^2), function(i) calc_mult_rowsum1(v = Index, u = Wtime22[, i], A = CondExp2, M = temp0a)) 
    }
    temp0c <- XZb2 *temp0a 
  }
  if (ncw > 0) temp5 <- lapply(1:(ncw), function(i) calc_mult_rowsum1(v = Index, u = Wtime2[, i], A = CondExp2, M = temp0c)) 

  Integral2 <- Integral[nk != 0,]
  post1 <- sum((temp1 * Integral2) %*% wGQ)
  post2 <- sum((temp2 * Integral2) %*% wGQ)
  if (ncw > 0) {
    post3 <- unlist(lapply(temp3, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw #
    post4 <- unlist(lapply(temp4, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw^2 #
    post5 <- unlist(lapply(temp5, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncw # 
    phiScore <- colSums(d * Wtime) - post3 # vector of length ncw #
  }
  
  if (model == 2) {
    alphaScore <- sum(d * rowSums(Ztime * post.bi)) - post1
  } else {
    alphaScore <- sum(d * Xtime %*% beta.old) + sum(d * rowSums(Ztime * post.bi)) - post1
  }
  if (ncw > 0) {
    pa.score <- c(phiScore, alphaScore)
    pa.info <- matrix(0, (ncw + 1), (ncw + 1)) # (ncw+1)*(ncw+1) matrix #
    pa.info[1:ncw, 1:ncw] <- -post4
    pa.info[(ncw + 1), 1:ncw] <- -post5
    pa.info[1:ncw, (ncw + 1)] <- -post5
    pa.info[(ncw + 1), (ncw + 1)] <- -post2
    
    #=============== Update phi and alpha ===============#
    pa.old <- c(phi.old, alpha.old) # vector of length (ncw+1) #
    paSVD <- svd(pa.info)
    pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
    pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncw+1) #
    phi.new <- pa.new[1:ncw]
    alpha.new <- pa.new[ncw + 1]
  } else {
    alpha.new <- alpha.old - alphaScore / (-post2)
    phi.new <- phi.old
  }

  Wtime2_phi.new <- if (ncw > 0) Wtime2 %*% phi.new else rep(0, M)
  if(model == 1) {
    #========== calculate the score of beta ==========#
    newZtime2.b = alpha.new * Ztime2.b
    eta.s.n1 <- as.vector(Wtime2_phi.new + alpha.new * Xtime2 %*% beta.old) + newZtime2.b # M*GQ matrix # 
    calc_expM2(eta.s.n1)
    temp0e <- alpha.new * eta.s.n1 * lamb.old[Index1]; 
    temp6 <- lapply(1:(ncx), function(i) calc_mult_rowsum1(v = Index, u = Xtime2[, i], A = CondExp2, M= temp0e))
    temp0d <- alpha.new * temp0e
    temp7 <- lapply(1:(ncx^2), function(i) calc_mult_rowsum1(v = Index, u = Xtime22[, i], A = CondExp2, M= temp0d))
    
    post6 <- unlist(lapply(temp6, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncx #
    post7 <- unlist(lapply(temp7, function(x) sum((x * Integral2) %*% wGQ))) # vector of length ncx^2 #
    post8 <- as.vector(Y - X %*% beta.old - rowSums(Z * post.bi[ID, ])) * X # N*ncx matrix #
    
    betaScore <- colSums(post8) / Ysigma2.new + alpha.new * colSums(d * Xtime) - post6
    betaInfo <- - X2.sum / Ysigma2.new - post7
    
    #========== Update beta ==========#
    betaSVD <- svd(betaInfo)
    betaInfo.inv <- betaSVD$v %*% diag(1 / betaSVD$d) %*% t(betaSVD$u)
    beta.new <- as.vector(beta.old - betaInfo.inv %*% betaScore) # vector of length ncx #
  }

  #========== Calculate the new lambda with new parameters ==========# 
  if ( model ==  2) { 
    eta.sn <- as.vector(Wtime2_phi.new) + alpha.new * Ztime2.b # M*GQ matrix # 
  } else {
    eta.sn <- as.vector(Wtime2_phi.new + alpha.new * Xtime2 %*% beta.new) + newZtime2.b  # M*GQ matrix # 
  }
  calc_expM2(eta.sn);
  calc_M1_M2_M3_Hadamard(eta.sn, CondExp ,  Integral, as.integer(Index - 1))
  tempLamb <- calc_M_v(v = wGQ, M = eta.sn)
  postLamb <- calc_tapply_vect_sum(v1 = tempLamb, v2 = as.integer(Index1 - 1)); ## Check this!

  lamb.new <- Index2 / postLamb
  
  result <- list(beta = beta.new, Ysigma = sqrt(Ysigma2.new), Bsigma = Bsigma.new, phi = phi.new, 
                 alpha = alpha.new, lamb = lamb.new, lgLik = lgLik, est.bi = post.bi)
  return(result)
}
