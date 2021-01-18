gpd.sandwich <-
    # Compute the filling in the Huber sandwich estimator of the covariance of gpd model parameters, by using the observed score vectors 
function(o){
    if (!inherits(o, "evmOpt")){ stop("object must be of class 'evmOpt'") }
    
    x <- o$data$D$phi; z <- o$data$D$xi
    ns <- ncol(x); nk <- ncol(z)
    phi <- coef(o)[1:ns]
    xi <- coef(o)[(ns+1):(ns + nk)]

    phi.i <- colSums(phi * t(x))
    xi.i <- colSums(xi * t(z))
    w.i <- (o$data$y - o$threshold) / exp(phi.i)

    if (any(xi.i < -.50)){ message("Fitted values of xi < -0.5") }

    # First derivatives of log-lik wrt coefficients of linear predictors

    dli.dphi <- (1 + 1/xi.i) * xi.i * w.i / (1 + xi.i*w.i) - 1
    dli.dxi <- 1/xi.i^2 * log(1 + xi.i*w.i)  - (1 + 1/xi.i)*w.i/(1 + xi.i*w.i)

    # Matrix Sc of score vectors, one row for each excess. 
    # First ns columns correspond to phi parameters, following nk columns correspond to xi

    nd <- nrow(x) # number of excesses
    Sc <- matrix(0, nrow=nd,ncol=ns+nk)
    for (s in 1:ns){
      Sc[,s] <- x[,s] * dli.dphi
    }
    for (k in 1:nk){
      Sc[,ns + k] <- z[,k] * dli.dxi
    }
    
    # now calculate observed covariance from observed scores for each data point and summing over the data:

    Cov.L1 <- matrix(0,nrow=ns+nk,ncol=ns+nk)
    for(u in 1:(ns+nk)){
        for(v in 1:(ns+nk)){
            Cov.L1[u,v] <- sum(Sc[,u] * Sc[,v])
        }
    }
    Cov.L1
}

