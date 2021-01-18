# Parameters


s2Params <- function(lambda1,
                       lambda2 = 0,
                       gamma1 = 0,
                       gamma2 = 0,
                       gamma3 = 0){
  s2P = c(
    lambda1 = lambda1,
    lambda2 = lambda2,
    gamma1 = gamma1,
    gamma2 = gamma2,
    gamma3 = gamma3
  )
  class(s2P) = "s2Params"
  return(s2P)
}

s2Fista <- function(MAX_ITER_INNER = 5000, TOL = 1e-7, t0 = 2, step = 0.1, use_warmstart = FALSE){
  s2F = list(
    MAX_ITER_INNER = MAX_ITER_INNER,
    TOL = TOL,
    t0 = t0,
    step = step,
    use_warmstart = use_warmstart
  )
  class(s2F) = "s2Fista"
  return(s2F)
}

print.s2Fista <- function(x, ...){
  print(unlist(x))
}

