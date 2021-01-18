fsfused <- function(y, s = 10, T, K.max=5){
  fit <- sl0fused_c(y, T0 = s, T02 = T, K.max)
  beta <- drop(fit$beta)
  d <- drop(fit$d)
  z <- drop(fit$z)
  u <- drop(fit$u)

  return(list(y = y, beta = beta, v = z))
}

