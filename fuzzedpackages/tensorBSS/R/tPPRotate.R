tPPRotate <- function(U0, xm, nl_int, eps, maxiter){
  
  delta <- eps + 1
  iter <- 0
  
  p <- nrow(U0)
  q <- ncol(xm)
  
  # The scale and the reference value
  scale_tau <- sum(diag(mModeCovariance(xm, 1, center = FALSE)))/p
  # Uses now only the x^2 objective function
  if(nl_int == 1){
    h2 <- scale_tau^2*q*(2 + q)
  }
  if(nl_int == 2){
    h2 <- scale_tau^(3/2)*(2^(3/2))*gamma(3/2 + q/2)/gamma(q/2)
  }
  if(nl_int == 3){
    f <- function(x) log(cosh(sqrt(x)))*(1/scale_tau)*dchisq(x/scale_tau, q)
    h2 <- integrate(Vectorize(f), 0, Inf)$value
  }
  
  while(delta > eps){
    Tmat <- matrix(0, p, p)
    for(i in 1:p){
      h1 <- computeh(U0[, i], xm, nl_int)
      mvec <- computeT(U0[, i], xm, nl_int)
      bk <- computeb(U0[, i], xm, nl_int)
      dk <- computed(U0[, i], xm, nl_int)
      Tmat[, i] <- (h1 - h2)*(mvec - scale_tau*(2*bk + q*dk)*U0[, i])
    }
    U1 <- Tmat%*%symmetricPower(crossprod(Tmat), -0.5)
    
    # Criterion evaluation
    J <- diag(sign(diag(t(U0)%*%U1)))
    
    delta <- sqrt(sum((U1 - U0%*%J)^2))
    U0 <- U1%*%J
    
    # if(verbose){
    #   print(delta)
    #   print(U1)
    # }
    
    iter <- iter + 1
    if(iter > maxiter){
      stop("Too many iterations")
    } 
  }
  return(list(U0, iter))
}