
#=============== Initial Value Calculation for Transformation Model II ===============#

InitValTMGeneric <- function (beta, model, n, X, Z, bBLUP, ID, Xtime, Ztime, Xtime2, Ztime2, Indcs, start, event, stop, W , ncw, Wtime2, rho, nk, Wtime22, d, Wtime, cvals) {

  Index = Indcs$Index  
  Index0 = Indcs$Index0
  Index1 = Indcs$Index1
  Index2 = Indcs$Index2

  tol.P = cvals$tol.P;
  iter = cvals$max.iter;
  
  M <- nrow(Xtime2)
   
  if(model == 1) { 
    fixedOrRand <- as.vector(X %*% beta) + rowSums(Z * bBLUP[ID, ]) # vector of length N #
    fixedOrRand.time <- as.vector(Xtime %*% beta) + rowSums(Ztime * bBLUP) # vector of length n #
    fixedOrRand.time2 <- as.vector(Xtime2 %*% beta) + rowSums(Ztime2 * bBLUP[Index, ]) # vector of length M #
  } else if(model == 2) {
    fixedOrRand <- rowSums(Z * bBLUP[ID, ]) # vector of length N #
    fixedOrRand.time <- rowSums(Ztime * bBLUP) # vector of length n #
    fixedOrRand.time2 <- rowSums(Ztime2 * bBLUP[Index, ]) # vector of length M # 
  } else {
    stop("Invalid model type")
  }
   
  #========== first fit the Cox model ==========#
  data.init <- data.frame(start = start, stop = stop, event = event, W = W, fixedOrRand = fixedOrRand)
  fit <- if (ncw > 0) coxph(Surv(start, stop, event) ~ W + fixedOrRand, data = data.init) else coxph(Surv(start, stop, event) ~ fixedOrRand, data = data.init)
  phi.old <- if (ncw > 0) fit$coefficients[1:ncw] else numeric(0)
  alpha.old <- fit$coefficients[ncw + 1]
  Wtime2_phi.old <- if (ncw > 0) Wtime2 %*% phi.old else rep(0, M)
  temp <- as.vector(exp(Wtime2_phi.old + alpha.old * fixedOrRand.time2)) # M*1 vector #
  lamb.old <- Index2 / calc_tapply_vect_sum(  v1 = temp, v2 = as.integer(Index1 - 1)) # vector of length n_u #
  
  if (rho == 0) {
    phi.new <- phi.old
    alpha.new <- alpha.old
    lamb.new <- lamb.old
  } else {
    for (it in 1:iter) {
      exp.es <- exp(as.vector(Wtime2_phi.old + alpha.old * fixedOrRand.time2))
      const <- rep(0, n)
      temp0a <- exp.es * lamb.old[Index1];
      const[nk != 0] <- calc_tapply_vect_sum(  v1 = temp0a, v2 = as.integer(Index - 1)) # vector of length n #
      
      CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|Oi), vector of length n #
      CondExp2 <- CondExp[nk != 0]
      
      temp0b <- fixedOrRand.time2 * temp0a;
      
      temp1 <- sum(CondExp2 *calc_tapply_vect_sum(  v1 = temp0b, v2 = as.integer(Index - 1)));
      temp2 <- sum(CondExp2 *calc_tapply_vect_sum(  v1 = fixedOrRand.time2 * temp0b, v2 = as.integer(Index - 1)));
      if (ncw > 0) {
        temp3 <- lapply(1:ncw, function(i) CondExp2 * calc_tapply_vect_sum(  v1 = Wtime2[, i] * temp0a, v2 =  as.integer(Index - 1)))
        temp3 <- sapply(temp3, sum) # vector of length ncw #
        temp4 <-  lapply(1:ncw^2, function(i) CondExp2 * calc_tapply_vect_sum(  v1 = Wtime22[, i] * temp0a, v2 =  as.integer(Index - 1)))
        temp4 <- sapply(temp4, sum) # vector of length ncw^2 #
        temp5 <- lapply(1:ncw, function(i) CondExp2 * calc_tapply_vect_sum(  v1 = Wtime2[, i] * temp0b, v2 =  as.integer(Index - 1)))
        temp5 <- sapply(temp5, sum) # vector of length ncw #  
        phiScore <- colSums(d * Wtime) - temp3 # vector of length ncw #
      }
      alphaScore <- sum(d * fixedOrRand.time) - temp1
      
      if (ncw > 0) {
        pa.score <- c(phiScore, alphaScore)
        pa.info <- matrix(0, (ncw + 1), (ncw + 1)) # (ncw+1)*(ncw+1) matrix #
        pa.info[1:ncw, 1:ncw] <- - temp4
        pa.info[(ncw + 1), 1:ncw] <- - temp5
        pa.info[1:ncw, (ncw + 1)] <- - temp5
        pa.info[(ncw + 1), (ncw + 1)] <- - temp2
        
        #=============== Update phi and alpha ===============#
        pa.old <- c(phi.old, alpha.old) # vector of length (ncw+1) #
        paSVD <- svd(pa.info)
        pa.info.inv <- paSVD$v %*% diag(1/paSVD$d) %*% t(paSVD$u)
        pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncw+1) #
        phi.new <- pa.new[1:ncw]
        alpha.new <- pa.new[ncw + 1]
      } else {
        alpha.new <- alpha.old - alphaScore / (-temp2)
        phi.new <- phi.old
        pa.new <- alpha.new
        pa.old <- alpha.old
      }
      
      Wtime2_phi.new <- if (ncw > 0) Wtime2 %*% phi.new else rep(0, M)
      #========== Calculate the new lambda with new parameters ==========#
      exp.esn <- exp(as.vector(Wtime2_phi.new + alpha.new * fixedOrRand.time2))
      tempLamb <- calc_tapply_vect_sum(  v1 = CondExp[Index] * exp.esn, v2 = as.integer(Index1 - 1))
      lamb.new <- Index2 / tempLamb
      
      #========== Check Convergence ==========#
      err <- max(abs(pa.new - pa.old) / (abs(pa.old) + tol.P)) 
      if (err <= tol.P) break
      else {
        phi.old <- phi.new
        alpha.old <- alpha.new
        lamb.old <- lamb.new
        Wtime2_phi.old <- if (ncw > 0) Wtime2 %*% phi.old else rep(0, M)
      }
    }
  }
  
  result <- list(phi = phi.new, alpha = alpha.new, lamb = lamb.new)
  return(result)
}
