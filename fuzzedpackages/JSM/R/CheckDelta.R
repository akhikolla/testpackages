
#========== Check whether the delta is too large when forward difference is used ==========#

CheckDeltaFD <- function (theta, ncz, delta) {
  
  Ysigma <- theta$Ysigma
  check1 <- (Ysigma - 2 * delta > 0)
  
  Bsigma <- theta$Bsigma
  if (is.matrix(Bsigma)) {
    VecB <- Bsigma[lower.tri(Bsigma, diag=TRUE)]
    check2 <- rep(FALSE, length(VecB))
    for (k in 1:length(VecB)) {
      temp.VecB <- VecB
      temp.VecB[k] <- VecB[k] + delta
      temp.B <- matrix(0, ncz, ncz)
      temp.B[lower.tri(temp.B, diag=TRUE)] <- temp.VecB
      temp.B <- temp.B + t(temp.B) - diag(diag(temp.B)) # ncz*ncz matrix #
      if(all(eigen(temp.B)$values > 0)) check2[k] <- TRUE
    }
  }else{
    check2 <- TRUE
  }
  return(all(c(check1, check2))) 
}


#========== Check whether the delta is too large when Richardson Extrapolation is used ==========#

CheckDeltaRE <- function (theta, ncz, delta) {
  
  Ysigma <- theta$Ysigma
  check1 <- (Ysigma - 2 * delta > 0)
  
  Bsigma <- theta$Bsigma
  if (is.matrix(Bsigma)) {
    VecB <- Bsigma[lower.tri(Bsigma, diag=TRUE)]
    check21 <- check22 <- rep(FALSE, length(VecB))
    for (k in 1:length(VecB)) {
      temp.VecB1 <- temp.VecB2 <- VecB
      temp.VecB1[k] <- VecB[k] - 2 * delta
      temp.VecB2[k] <- VecB[k] + 2 * delta
      temp.B1 <- temp.B2 <- matrix(0, ncz, ncz)
      temp.B1[lower.tri(temp.B1, diag=TRUE)] <- temp.VecB1
      temp.B2[lower.tri(temp.B2, diag=TRUE)] <- temp.VecB2
      temp.B1 <- temp.B1 + t(temp.B1) - diag(diag(temp.B1)) # ncz*ncz matrix #
      temp.B2 <- temp.B2 + t(temp.B2) - diag(diag(temp.B2)) # ncz*ncz matrix #
      if(all(eigen(temp.B1)$values > 0)) check21[k] <- TRUE
      if(all(eigen(temp.B2)$values > 0)) check22[k] <- TRUE
    }
  }else{
    check21 <- check22 <- (Bsigma - 2 * delta > 0)
  }
  return(all(c(check1, check21, check22))) 
}

