# this function samples priors for the estimates and the number of indicators

priorSamp <- function(p, estimates, n.samp = 2e3){
  if ("alpha" %in% estimates || "lambda2" %in% estimates || "lambda4" %in% estimates || "lambda6" %in% estimates ||
      "glb" %in% estimates){
    v0 <- p
    k0 <- 1e-10
    t <- diag(p)
    T0 <- solve(t/k0)
    m <- array(0, c(n.samp, p, p))
    for (i in 1:n.samp){
      m[i, , ] <- LaplacesDemon::rinvwishart(v0, T0)
    }
  }
  out <- list()
  if ("alpha" %in% estimates){
    priora <- apply(m, MARGIN = 1, applyalpha)
    out$prioralpha <- quantiles(priora[priora >= 0])
  }
  if ("lambda2" %in% estimates){
    priorlambda2 <- apply(m, MARGIN = 1, applylambda2)
    out$priorlambda2 <- quantiles(priorlambda2[priorlambda2 >= 0])
  }
  if ("lambda4" %in% estimates){
    priorlambda4 <- apply(m, MARGIN = 1, applylambda4_nocpp)
    out$priorlambda4 <- quantiles(priorlambda4[priorlambda4 >= 0])
  }
  if ("lambda6" %in% estimates){
    priorlambda6 <- apply(m, MARGIN = 1, applylambda6)
    out$priorlambda6 <- quantiles(priorlambda6[priorlambda6 >= 0])
  }
  if ("glb" %in% estimates){
    # control <- Rcsdp::csdp.control(printlevel = 0)
    # write.control.file(control)
    priorglb <- glbOnArray_custom(m)
    out$priorglb <- quantiles(priorglb[priorglb >= 0])
    # unlink("param.csdp")
  }
  if ("omega" %in% estimates){
    H0 <- 1 # prior multiplier matrix for lambdas variance
    l0k <- rep(0, p) # prior lambdas
    a0k <- 1 # prior gamma function for psis
    b0k <- 2 # prior gamma for psi
    prioromega <- numeric(n.samp)
    for (i in 1:n.samp){
      invpsi <- rgamma(p, a0k, b0k)
      invPsi <- diag(invpsi)
      psi <- 1/invpsi
      lambda <- rnorm(p, l0k, sqrt(psi * H0))
      prioromega[i] <- omegaBasic(lambda, psi)
    }
    out$prioromega <- quantiles(prioromega[prioromega >= 0])
  }

  return(out)

}


