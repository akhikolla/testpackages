
#=============== Function to Calculate the Likelihood Value for Model II ===============#
#LHGeneric <- function (theta){
LHGeneric <- function (theta, n, Z.st, Y.st, X.st, b, Ztime, nk, Index0, model, Wtime, Xtime, Wtime2, Xtime2, GQ, Index, Index1, rho, d, wGQ, ncz, Ztime2.st, ncw) {
  
  beta <- theta$beta
  Ysigma2 <- (theta$Ysigma) ^ 2
  Bsigma <- theta$Bsigma
  phi <- theta$phi
  alpha <- theta$alpha
  lamb <- theta$lamb
  
  M <- nrow(Xtime2)

  VY <- lapply(1:n, function(i) calc_VY(M = Z.st[[i]], A = Bsigma, b = Ysigma2 ))  
  VB <-  lapply(1:n, function(i) calc_VB( Bsigma ,M2 =  Z.st[[i]], M3 = VY[[i]])) 
  muB <- lapply(1:n, function(i) calc_muB( BSold=Bsigma , Zst=Z.st[[i]], Yst=Y.st[[i]], betaold=beta ,VY= VY[[i]], Xst=X.st[[i]]))
  bi.st <- lapply(1:n, function(i) calc_bi_st(v0=muB[[i]], b ,M = VB[[i]]) ) 

  # each element is ncz*GQ matrix #
  bi <- do.call(rbind, bi.st)
  Ztime.b <- do.call(rbind, lapply(1:n, function(i) Ztime[i, ] %*% bi.st[[i]])) # n*GQ matrix # 
  Ztime2.b <-fast_lapply_length(Ztime2.st, bi.st, (1:n)[nk !=      0] - 1)# M*GQ matrix # 

  log.lamb <- log(lamb[Index0])
  log.lamb[is.na(log.lamb)] <- 0  
  
  Wtime_phi <- if (ncw > 0) Wtime %*% phi else rep(0, n)
  Wtime2_phi <- if (ncw > 0) Wtime2 %*% phi else rep(0, M)
  
  if (model == 2) { 
    log.density1 <- log.lamb + as.vector(Wtime_phi) + alpha * Ztime.b # n*GQ matrix #
  } else if(model ==1){ 
    log.density1 <- log.lamb + as.vector(Wtime_phi + alpha * Xtime %*% beta) + alpha * Ztime.b # n*GQ matrix # 
  } else {
    stop("Invalid model type")
  }

  calc_v_a( Ztime2.b, alpha); # Ztime2.b gets altered

  if ( model == 2) {
    exp.es<- as.numeric(Wtime2_phi) + Ztime2.b  
  } else {
    exp.es<- as.numeric(Wtime2_phi + alpha * Xtime2 %*% beta) + Ztime2.b  
  }
  calc_expM2(exp.es) 

  const <- matrix(0, n, GQ) # n*GQ matrix # 
  const[nk != 0, ] <- calc_rowsum( (Index),  exp.es * lamb[Index1])

  log.density2 <- - log(1 + rho * const) # n*GQ matrix #  
  log.survival <- if(rho > 0) log.density2 / rho else - const # n*GQ matrix #   
  f.surv <- exp(d * log.density1 + d * log.density2 + log.survival) # n*GQ matrix #
  deno <- as.vector(f.surv %*% wGQ) # vector of length n #
   
  f.long <- sapply(1:n, function(i) calc_MVND(Y.st[[i]], as.vector(X.st[[i]] %*% beta), VY[[i]]))

  # vector of length n #
  return(sum(log(f.long * deno / (pi ^ (ncz / 2)))))
}
