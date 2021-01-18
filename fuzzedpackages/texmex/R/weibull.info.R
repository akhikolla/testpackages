weibull.info <-
  # phi=log(beta) and gamma are both linear in their covariates. 
  # If penalization is used, the calculation accounts for this, but the resulting
  # estimates of variance will be too low and bias might dominate MSE
  function(o, method="observed"){
    if (!inherits(o, "evmOpt")){ stop("object must be of class 'evmOpt'") }
    if (method != "observed"){ stop("only 'observed' information is implemented") }
  
    x <- o$data$D$phi
    z <- o$data$D$gamma
    ns <- ncol(x)
    nk <- ncol(z)
    phi <- coef(o)[1:ns]
    gamma <- coef(o)[(ns+1):(ns + nk)]
    
    phi.i <- colSums(phi * t(x))
    gamma.i <- colSums(gamma * t(z))
    w.i <- (o$data$y - o$threshold)
    
    # Second derivatives of penalties
    p <- matrix(0, nrow=ns+nk, ncol=ns+nk)
    if (o$penalty %in% c("gaussian", "quadratic")){ # note if Lasso penalty used then 2nd deriv is zero hence no term for this
      Si <- solve(o$priorParameters[[2]])
      for (i in 1:(ns+nk)){
        for (j in 1:(ns + nk)){
          p[i,j] <- 2*Si[i,j]
        }
      }
    }
    
    # Second and mixed derivatives of log-lik wrt coefficients of linear predictors
    
    # linear predictors in beta  
#    d2li.dbeta2 <- gamma.i / beta.i^2 * (1 - (1 + gamma.i) * (w.i/beta.i)^gamma.i )
#    d2li.dbetadgamma <- 1/beta.i  * ( (w.i/beta.i)^gamma.i - 1) + gamma.i/beta.i * (w.i/beta.i)^gamma.i / w.i * log(w.i/beta.i)
#    d2li.dgamma2 <- - gamma.i^(-2) - (w.i/beta.i)^gamma.i * (log(w.i/beta.i))^2
    
    d2li.dphi2 <- -gamma.i  * (w.i^gamma.i ) * exp(-gamma.i*phi.i) * (1 + gamma.i)
    d2li.dphidgamma <- w.i^gamma.i * exp( - gamma.i * phi.i ) *  (gamma.i * (log(w.i) - phi.i) + 1 ) - 1
    d2li.dgamma2 <- - gamma.i^(-2) - (w.i/exp(phi.i))^gamma.i * (log(w.i) - phi.i)^2
    
    # Matrix has 4 blocks, 2 of which are transposes of each other. Need block for beta parameters,
    # block for gamma parameters and block for the cross of them.
    
    Ip <- matrix(0, ncol=ns, nrow=ns)
    for (u in 1:ns){
      for (v in 1:ns){
        Ip[u,v] <- -sum(x[,u] * x[,v] * d2li.dphi2)
      }
    }
    
    Ix <- matrix(0, ncol=nk, nrow=nk)
    for (s in 1:nk){
      for (t in 1:nk){
        Ix[s,t] <- -sum(z[,s] * z[,t] * d2li.dgamma2)
      }
    }
    
    Ipx <- matrix(0, ncol=nk, nrow=ns)
    for (u in 1:ns){
      for (s in 1:nk){
        Ipx[u,s] <- -sum(z[,s] * x[,u] * d2li.dphidgamma )
      }
    }

        i <- rbind( cbind(Ip, Ipx), cbind(t(Ipx), Ix))
    
    # return observed Information matrix.   Note that an estimate of the covariance matrix is given by the inverse of this matrix.
    i - p
  }
