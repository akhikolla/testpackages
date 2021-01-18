

### y      - numeric, dependent variable
### x      - matrix 1 x m of independent variables
### m      - number of independent variables
### theta  - matrix 1 x m of regression coefficients
### E      - state covariance matrix 
### V      - output variance
### lambda - forgetting factor
### kappa  - parameter in Exponentially Weighted Moving Averaging
### t      - time index
  
  
.kalman <- function(y,x,theta,E,V,lambda,kappa,t) 
  {
    yhat <- as.numeric(x %*% t(theta))
    e <- y - yhat
    R <- E / lambda
    Rf <- sqrtmat(R)
    Vt <- as.numeric((x %*% Rf) %*% t(x %*% Rf))
    Vu <- V + Vt
    E <- R - (R %*% t(x)) %*%  (x %*% R) / Vu
    pdens <- exp(-0.5 * e^2 / Vu) / sqrt(2 * pi * Vu)
    theta <- theta + t((R %*% t(x)) * e / Vu)
    if (!is.null(kappa))
      {
        V <- V * kappa + (1 - kappa) * e^2
      }
    else
      {
        temp <- ((t-1) * V + e^2 - Vt) / t
        if (temp > 0) { V <- temp }
      }
    
    out <- list(yhat,theta,E,V,pdens)
    names(out) <- c("y.hat","theta","E","V","pdens")
    return(out)
  }


.kalman2 <- function(y,x,theta,R,t,Rw,Vv) 
  {
    yhat <- as.numeric(x %*% t(theta))
    e <- y - yhat
    R <- Rw + R
    Rf <- sqrtmat(R)
    Vt <- as.numeric((x %*% Rf) %*% t(x %*% Rf))
    Vu <- Vv + Vt
    R <- R - (R %*% t(x)) %*%  (x %*% R) / Vu
    pdens <- exp(-0.5 * e^2 / Vu) / sqrt(2 * pi * Vu)
    theta <- theta + t((R %*% t(x)) * e / Vv)
    
    out <- list(yhat,theta,R,Vu,pdens)
    names(out) <- c("y.hat","theta","E","V","pdens")
    return(out)
  }
  
  