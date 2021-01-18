
ffused <- function(y, s = 10, K.max=5){

  fit <- l0fused_c(y, s, K.max)
  beta <- drop(fit$beta)
  z <- drop(fit$z)
  u <- drop(fit$u)
  return(list(y = y, beta = beta, v = z))
}
