
#=============== Initial Value Calculation for Model II with NMRE ===============#

InitValMultGeneric <- function (gamma, B.st, n, Y.st, ni, model, ID, Index, B, Btime, Btime2, start, stop, event, Z, ncz, Ztime2, Index2, Index1, rho, iter, nk, d, Ztime22, Ztime, tol.P) {
  
  BTg <- lapply(B.st, function(x) as.vector(x %*% gamma))
  G <- unlist(lapply(1:n, function(i)  as.vector(   tcrossprod(Y.st[[i]] -  BTg[[i]]))))
  tempCov1 <- unlist(lapply(BTg, function(x) tcrossprod(x) ))
  tempCov2 <- unlist(lapply(1:n, function(i) c(diag(1, ni[i]))))
  tempCov <- cbind(tempCov1, tempCov2)
  sigmas <- solve(t(tempCov) %*% tempCov) %*% t(tempCov) %*% G
  Bsigma2 <- sigmas[1]
  Ysigma2 <- sigmas[2]
  if(Bsigma2 < 0) Bsigma2 <- 0.1
  if(Ysigma2 < 0) Ysigma2 <- 0.1
  
  M <- length(Index)
  
  VY <- lapply(1:n, function(i) calc_VY( M = BTg[[i]], A = Bsigma2, b = Ysigma2))  
  bBLUP <-unlist(lapply(1:n, function(i) calc_muBMult(  Bsigma2,VY[[i]],BTg[[i]],Y.st[[i]] )+1 ))
  if (model == 2){
    fixedOrRand <- bBLUP[ID]
    fixedOrRand.time2 <- bBLUP[Index]
  } else if( model == 1) {
    fixedOrRand <- bBLUP[ID] * as.vector(B %*% gamma)
    fixedOrRand.time <- bBLUP * as.vector(Btime %*% gamma)
    fixedOrRand.time2 <- bBLUP[Index] * as.vector(Btime2 %*% gamma)
  } else {
    stop("Invalid model type")     
  }
  
  #========== first fit the Cox model ==========#
  data.init <- data.frame(start = start, stop = stop, event = event, Z = Z, fixedOrRand = fixedOrRand)
  fit <- if (ncz > 0) coxph(Surv(start, stop, event) ~ Z + fixedOrRand, data = data.init) else coxph(Surv(start, stop, event) ~ fixedOrRand, data = data.init)
  phi.old <- if (ncz > 0) fit$coefficients[1:ncz] else numeric(0)
  alpha.old <- fit$coefficients[ncz + 1]
  Ztime2_phi.old <- if (ncz > 0) Ztime2 %*% phi.old else rep(0, M)
  temp <- as.vector(exp(Ztime2_phi.old + alpha.old * fixedOrRand.time2)) # M*1 vector #
  lamb.old <- Index2 / calc_tapply_vect_sum( temp,  as.integer(Index1-1))

  if (rho == 0) {
    phi.new <- phi.old
    alpha.new <- alpha.old
    lamb.new <- lamb.old
  } else {
   for (it in 1:iter) { 
      exp.es <- exp(as.vector(Ztime2_phi.old + alpha.old * fixedOrRand.time2))
      temp0a <- exp.es * lamb.old[Index1];
      temp0b <- fixedOrRand.time2 * temp0a; 
      const <- rep(0, n)
      const[nk != 0] <- calc_tapply_vect_sum(v1 = temp0a, v2 = as.integer(Index - 1)) # vector of length n #
      CondExp <- (1 + d * rho) / (1 + rho * const) # conditional expectation E(xi|Oi), vector of length n #
      CondExp2 <- CondExp[nk != 0]
      
      if (ncz > 0) {
        temp1 <- lapply(1:ncz, function(i) CondExp2 * calc_tapply_vect_sum(v1 = Ztime2[, i] * temp0a, v2 = as.integer(Index - 1)))
        temp1 <- sapply(temp1, sum) # vector of length ncz #
        temp3 <- lapply(1:(ncz ^ 2), function(i) CondExp2 * calc_tapply_vect_sum(v1 = Ztime22[, i] * temp0a, v2 = as.integer(Index - 1)))
        temp3 <- sapply(temp3, sum) # vector of length ncz^2 #
        temp5 <- lapply(1:ncz, function(i) CondExp2 * calc_tapply_vect_sum(v1 = Ztime2[, i] * temp0b, v2 = as.integer(Index - 1)))
        temp5 <- sapply(temp5, sum) # vector of length ncz #
        phiScore <- colSums(d * Ztime) - temp1 # vector of length ncz #
      }
      temp2 <- sum(CondExp2 * calc_tapply_vect_sum(v1 = temp0b, v2 = as.integer(Index - 1)))
      temp4 <- sum(CondExp2 * calc_tapply_vect_sum(v1 = fixedOrRand.time2 * temp0b, v2 = as.integer(Index - 1)))

      if (model==2) {
        alphaScore <- sum(d * bBLUP) - temp2
      } else {
        alphaScore <- sum(d * fixedOrRand.time) - temp2
      }
      
      if(ncz > 0) {
        pa.score <- c(phiScore, alphaScore)
        pa.info <- matrix(0, (ncz + 1), (ncz + 1)) # (ncz+1)*(ncz+1) matrix #
        pa.info[1:ncz, 1:ncz] <- - temp3
        pa.info[(ncz + 1), (ncz + 1)] <- - temp4
        pa.info[(ncz + 1), 1:ncz] <- - temp5
        pa.info[1:ncz, (ncz + 1)] <- - temp5
        
        #=============== Update phi and alpha ===============#
        pa.old <- c(phi.old, alpha.old) # vector of length (ncz+1) #
        paSVD <- svd(pa.info)
        pa.info.inv <- paSVD$v %*% diag(1 / paSVD$d) %*% t(paSVD$u)
        pa.new <- pa.old - pa.info.inv %*% pa.score # vector of length (ncz+1) #
        phi.new <- pa.new[1 : ncz]
        alpha.new <- pa.new[ncz + 1]
      } else {
        alpha.new <- alpha.old - alphaScore / (-temp4)
        phi.new <- phi.old
        pa.new <- alpha.new
        pa.old <- alpha.old
      }
      
      Ztime2_phi.new <- if (ncz > 0) Ztime2 %*% phi.new else rep(0, M)
      #========== Calculate the new lambda with new parameters ==========#
      exp.esn <- exp(as.vector(Ztime2_phi.new + alpha.new * fixedOrRand.time2))
      tempLamb <- calc_tapply_vect_sum(v1 = CondExp[Index] * exp.esn, v2 = as.integer(Index1 - 1))
      lamb.new <- Index2 / tempLamb
      
      #========== Check Convergence ==========#
      err <- max(abs(pa.new - pa.old) / (abs(pa.old) + tol.P))
      if (err <= tol.P) break
      else {
        phi.old <- phi.new
        alpha.old <- alpha.new
        lamb.old <- lamb.new
        Ztime2_phi.old <- if (ncz > 0) Ztime2 %*% phi.old else rep(0, M) # Added by Pantelis
      }
    }
  }
  
  result <- list(phi = phi.new, alpha = alpha.new, lamb = lamb.new, Ysigma = sqrt(Ysigma2), Bsigma = sqrt(Bsigma2))
  return(result)
}
