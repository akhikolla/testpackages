
#=============== The DQ Function for Model II ===============#
#=============== Transformation model is fitted for the survival part ===============#

DQfuncGeneric <- function (model, ptheta, theta, n, Z.st, Y.st, X.st, Ztime, nk, Wtime, Ztime.b, Wtime2, Xtime, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncx, ncw, p, ncz, ncz2, b, Ztime2.st, Index0, X, Y, ID, N, Index2, Z) { # ptheta means "theta prime"
  ## This might have significant (e-14) numerical difference with the original code. Something odd is happening with VB/muB.
  pbeta <- ptheta$beta
  beta <- theta$beta
  pYsigma2 <- (ptheta$Ysigma) ^ 2
  Ysigma2 <- (theta$Ysigma) ^ 2
  pBsigma <- ptheta$Bsigma
  Bsigma <- theta$Bsigma
  pphi <- ptheta$phi
  phi <- theta$phi
  palpha <- ptheta$alpha
  alpha <- theta$alpha
  plamb <- ptheta$lamb
  lamb <- theta$lamb
  
  M <- nrow(Xtime2)
  
  VY <- lapply(1:n, function(i) calc_VY(M = Z.st[[i]], A = Bsigma, b = Ysigma2 ))  
  VB <-  lapply(1:n, function(i) calc_VB( Bsigma ,M2 =  Z.st[[i]], M3 = VY[[i]])) 
  muB <- lapply(1:n, function(i) calc_muB( BSold=Bsigma , Zst=Z.st[[i]], Yst=Y.st[[i]], betaold=beta ,VY= VY[[i]], Xst=X.st[[i]]))
  bi.st <- lapply(1:n, function(i) calc_bi_st(v0=muB[[i]], b ,M = VB[[i]]) ) 

  bi <- do.call(rbind, bi.st) # (n*ncz)*GQ matrix #
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix # 

  Ztime2.b <-fast_lapply_length(Ztime2.st, bi.st, (1:n)[nk != 0] - 1) # M*GQ matrix # 
  
  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0
  
  Wtime_phi <- if (ncw > 0) Wtime %*% phi else rep(0, n)
  Wtime2_phi <- if (ncw > 0) Wtime2 %*% phi else rep(0, M)
  Wtime2_pphi <- if (ncw > 0) Wtime2 %*% pphi else rep(0, M)

  if (model == 2) {
    log.density1 <- log.lamb + as.vector(Wtime_phi) + alpha * Ztime.b # n*GQ matrix #
    exp.es <- as.numeric(Wtime2_phi) + alpha * Ztime2.b  
  } else if (model == 1) {
    log.density1 <- log.lamb + as.vector(Wtime_phi + alpha * Xtime %*% beta) + alpha * Ztime.b # n*GQ matrix # 
    exp.es <- as.vector(Wtime2_phi + alpha * Xtime2 %*% beta) + alpha*Ztime2.b  
  } else {
    stop("Invalid model type")
  } 
  calc_expM2(exp.es)

  const <- matrix(0, n, GQ) # n*GQ matrix #

  const[nk != 0, ] <- calc_rowsum_mult((Index), lamb[Index1], exp.es)
  log.density2 <- - log(1 + rho * const) # n*GQ matrix # 
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*GQ matrix #
  
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
  Integral <- f.surv / deno # n*GQ matrix #
  CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|bi,Oi), n*GQ matrix #
  
  len <- ncx + ncw + p + 2
  Q <- rep(0, len)
  
  post.bi <- Integral %*% (t(bi) * wGQ) # n*(n*ncz) matrix #
  post.bi <- if(ncz > 1) t(sapply(1:n, function(i) post.bi[i, ((i - 1) * ncz + 1):(i * ncz)])) else 
             matrix(diag(post.bi), nrow = n) # n*ncz matrix #
  
  if (ncz > 1) {    
    tempB <-  fast_rbind_lapply_outerprod( bi.st )  # (n*ncz^2)*GQ matrix #     
  } else {
    tempB <- bi ^ 2
  }
  post.bi2 <- Integral %*% (t(tempB) * wGQ) # n*(n*ncz^2) matrix #
  post.bi2 <- if(ncz > 1) t(sapply(1:n, function(i) post.bi2[i, ((i - 1) * ncz2 + 1):(i * ncz2)])) else 
              matrix(diag(post.bi2), nrow = n) # n*(ncz^2) matrix #
  
  pYmu <- as.vector(X %*% pbeta) + do.call(rbind, lapply(1:n, function(i) Z.st[[i]] %*% bi.st[[i]])) 
  post.resid <- ((Y - pYmu) ^ 2 * Integral[ID, ]) %*% wGQ # vector of length N #
  Q[ncx + ncw + 2] <- - N / sqrt(pYsigma2) + sum(post.resid) / (pYsigma2 ^ (3 / 2))
  
  pBsigmaInv = solve(pBsigma);
  tempB <- - n * pBsigmaInv / 2 + pBsigmaInv %*% matrix(colSums(post.bi2), ncz, ncz) %*% pBsigmaInv / 2

  ind <- Indexing(ncz)
  Q[(ncx + ncw + 3):len] <- as.vector(tapply(c(tempB), ind, sum))

  if (model == 2) {
    exp.esp <-as.numeric(Wtime2_pphi) + alpha *Ztime2.b # M*GQ matrix #
  } else if(model == 1) {
    exp.esp <- as.vector(Wtime2_pphi + palpha * Xtime2 %*% pbeta) +  alpha *Ztime2.b # M*GQ matrix #
    XZb2 <- as.vector(Xtime2 %*% pbeta) + Ztime2.b # M*GQ matrix # 
  }  
  calc_expM2(exp.esp)

  temp0 <- exp.esp; temp0[1] = temp0[1] + 0 # "touch the variable"
  calc_M1_M2_M3_Hadamard(temp0, CondExp, Integral, as.integer(Index - 1))
  temp1 <- calc_M_v(v = wGQ, M = temp0)
  if (model == 2) {
    calc_M1_M2_Hadamard(temp0, Ztime2.b) }
  else { 
    calc_M1_M2_Hadamard(temp0, XZb2)
  }
  temp2 <- calc_M_v(v = wGQ, M = temp0)
  post1 <- calc_tapply_vect_sum( temp1, as.integer(Index1 - 1));
  post2 <- calc_tapply_vect_sum( temp2, as.integer(Index1 - 1));

  if (ncw > 0) {
    temp3 <- Wtime2 * temp1 # M*ncw matrix #
    post3 <- as.matrix(apply(temp3, 2, function(x) calc_tapply_vect_sum( x, as.integer(Index1-1))))
    Q[(ncx + 1):(ncx + ncw)] <- colSums(d * Wtime) - colSums(Index2 * post3 / post1) # vector of length ncw #
  }  

  if (model == 2) {    
    Q[ncx + ncw + 1] <- sum(d * rowSums(Ztime * post.bi)) - sum(Index2 * post2 / post1)
    pResid <- as.vector(Y - X %*% pbeta) - rowSums(Z * post.bi[ID, ]) # vector of length N #
    Q[1:ncx] <- colSums(X * pResid) / pYsigma2 # vector of length ncx # 
  } else { 
    Q[ncx + ncw + 1] <- sum(d * Xtime %*% pbeta) + sum(d * rowSums(Ztime * post.bi)) - sum(Index2 * post2 / post1)
    temp4 <- Xtime2 * temp1 # M*ncx matrix #
    post4 <- palpha * as.matrix(apply(temp4, 2, function(x) calc_tapply_vect_sum( x, as.integer(Index1-1)) )) # n_u*ncx matrix #
    pResid <- as.vector(Y - X %*% pbeta) - rowSums(Z * post.bi[ID, ]) # vector of length N #
    Q[1:ncx] <- colSums(X * pResid) / pYsigma2 + palpha * colSums(d * Xtime) - colSums(Index2 * post4 / post1) 
  }  

return(Q)
}
