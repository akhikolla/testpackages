
ftf <- function(y, k = 1, s = 20, K.max=5){

  fit = l0tf_c(y, k, s, K.max)
  beta = drop(fit$beta)
  v = drop(fit$z)
  u = drop(fit$u)
  return(list(y = y, beta = beta, v = v))
}

