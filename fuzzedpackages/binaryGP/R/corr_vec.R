corr_vec <- function(x, X, rho, corr){
  n <- nrow(X)
  d <- ncol(X)
  if (corr$type == "exponential") {
      r <- exp(-(abs(X - as.matrix(rep(1, n)) %*% (x))^corr$power) %*%
                (10^rho))

  }else if (corr$type == "matern") {
      temp <- 10^rho
      temp <- matrix(temp, ncol = d, nrow = (length(X)/d), byrow = TRUE)
      temp <- 2 * sqrt(corr$nu) * abs(X - as.matrix(rep(1, n)) %*% (x)) * (temp)
      ID <- which(temp == 0)
      rd <- (1/(gamma(corr$nu) * 2^(corr$nu - 1))) * (temp^corr$nu) * besselK(temp, corr$nu)
      rd[ID] <- 1
      r <- matrix(apply(rd, 1, prod), ncol = 1)
  }
  return(t(r))
}
