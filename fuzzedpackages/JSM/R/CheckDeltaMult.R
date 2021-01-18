
#========== Check whether the delta is too large when Richardson Extrapolation is used ==========#
#========== Multiplicative Joint Modeling ==========#

CheckDeltaMult <- function (theta, delta) {
  
  Ysigma <- theta$Ysigma
  check1 <- (Ysigma - 2 * delta > 0)
  
  Bsigma <- theta$Bsigma
  check2 <- (Bsigma - 2 * delta > 0)

  return(all(c(check1, check2))) 
}
