
roll.reg <- function(y,x=NULL,window,c=NULL)
{

### estimates rolling regression model

### y - a numeric or a column matrix of a dependent variable

### x - a matrix of independent variables (drivers), different columns correspond to different variables

### window - a numeric, a size of a window for rolling

### c - a parameter indicating whether constant is included,
###     by default c=TRUE (constant is included),
###     it is not possible to set c=FALSE if x=NULL 

if (is.null(c)) { c <- TRUE }

if (is.null(x)) { c <- TRUE }

if (! is.null(x)) { x <- as.matrix(x) }

y.roll.ols <- vector()
aic <- vector()
aicc <- vector()
bic <- vector() 
mse <- vector() 
if (! is.null(x))
  {
    if (c==TRUE)
      {
        coeff <- matrix(NA,ncol=(ncol(x)+1),nrow=1)
      }
    else
      {
        coeff <- matrix(NA,ncol=ncol(x),nrow=1)
      }
  }
else
  {
    coeff <- matrix(NA,ncol=1,nrow=1)
  }
pval <- coeff

if (! is.null(x))
  {
    for (i in 1:(window-1))
      {
        if (i==1) 
          { 
            if (c==TRUE) 
              {
                m <- lm(y[1] ~ t(x[1,]))
                n.par <- ncol(x) + 2
              }
            else
              {
                m <- lm(y[1] ~ t(x[1,]) -1)
                n.par <- ncol(x) + 1
              }
            y.roll.ols[1] <- m$fitted.values[1]
            aic[1] <- AIC(m)
            aicc[1] <- aic[1] + (2*n.par*(n.par+1))/(1-n.par-1)
            bic[1] <- BIC(m)
            mse[1] <- (m$residuals)^2
            mm <- summary(m)
            if (all(is.finite(mm$coefficients[,1])) && all(is.finite(mm$coefficients[,4])))
              {
                coeff <- rbind(coeff,mm$coefficients[,1])
                pval <- rbind(pval,mm$coefficients[,4])
              }
            else
              {
                coeff <- rbind(coeff,rep.int(0,ncol(coeff)))
                pval <- rbind(pval,rep.int(1,ncol(pval)))
              }
          }
        else
          {
            if (c==TRUE) 
              {
                m <- lm(y[1:i] ~ x[1:i,])
                n.par <- ncol(x) + 2
              }
            else
              {
                m <- lm(y[1:i] ~ x[1:i,] -1)
                n.par <- ncol(x) + 1
              }
            y.roll.ols[i] <- m$fitted.values[i]
            aic[i] <- AIC(m)
            aicc[i] <- aic[i] + (2*n.par*(n.par+1))/(i-n.par-1)
            bic[i] <- BIC(m)
            mse[i] <- mean((m$residuals)^2)
            mm <- summary(m)
            if (all(is.finite(mm$coefficients[,1])) && all(is.finite(mm$coefficients[,4])))
              {
                coeff <- rbind(coeff,mm$coefficients[,1])
                pval <- rbind(pval,mm$coefficients[,4])
              }
            else
              {
                coeff <- rbind(coeff,rep.int(0,ncol(coeff)))
                pval <- rbind(pval,rep.int(1,ncol(pval)))
              }
          }
      }

    for (i in window:nrow(x))
      {
        if (c==TRUE) 
          {
            m <- lm(y[(i-window+1):i] ~ x[(i-window+1):i,])   
            n.par <- ncol(x) + 2           
          }
        else
          {
            m <- lm(y[(i-window+1):i] ~ x[(i-window+1):i,] -1)  
            n.par <- ncol(x) + 1            
          }
        y.roll.ols[i] <- m$fitted.values[window]
        aic[i] <- AIC(m)
        aicc[i] <- aic[i] + (2*n.par*(n.par+1))/(window-n.par-1)
        bic[i] <- BIC(m)
        mse[i] <- mean((m$residuals)^2)
        mm <- summary(m)
        if (all(is.finite(mm$coefficients[,1])) && all(is.finite(mm$coefficients[,4])))
          {
            coeff <- rbind(coeff,mm$coefficients[,1])
            pval <- rbind(pval,mm$coefficients[,4])
          }
        else
          {
            coeff <- rbind(coeff,rep.int(0,ncol(coeff)))
            pval <- rbind(pval,rep.int(1,ncol(pval)))
          }
      }
  }
else
  {
    for (i in 1:(window-1))
      {
        if (i==1) 
          { 
            m <- lm(y[1] ~ 1)
            y.roll.ols[1] <- m$fitted.values[1]
            aic[1] <- AIC(m)
            n.par <- 2
            aicc[1] <- aic[1] + (2*n.par*(n.par+1))/(1-n.par-1)
            bic[1] <- BIC(m)
            mse[1] <- (m$residuals)^2
            mm <- summary(m)
            if (all(is.finite(mm$coefficients[,1])) && all(is.finite(mm$coefficients[,4])))
              {
                coeff <- rbind(coeff,mm$coefficients[,1])
                pval <- rbind(pval,mm$coefficients[,4])
              }
            else
              {
                coeff <- rbind(coeff,rep.int(0,ncol(coeff)))
                pval <- rbind(pval,rep.int(1,ncol(pval)))
              }
          }
        else
          {
            m <- lm(y[1:i] ~ 1)
            y.roll.ols[i] <- m$fitted.values[i]
            aic[i] <- AIC(m)
            n.par <- 2
            aicc[i] <- aic[i] + (2*n.par*(n.par+1))/(i-n.par-1)
            bic[i] <- BIC(m)
            mse[i] <- mean((m$residuals)^2)
            mm <- summary(m)
            if (all(is.finite(mm$coefficients[,1])) && all(is.finite(mm$coefficients[,4])))
              {
                coeff <- rbind(coeff,mm$coefficients[,1])
                pval <- rbind(pval,mm$coefficients[,4])
              }
            else
              {
                coeff <- rbind(coeff,rep.int(0,ncol(coeff)))
                pval <- rbind(pval,rep.int(1,ncol(pval)))
              }
          }
      }

    for (i in window:length(as.vector(y)))
      {
        m <- lm(y[(i-window+1):i] ~ 1)
        y.roll.ols[i] <- m$fitted.values[window]
        aic[i] <- AIC(m)
        n.par <- 2
        aicc[i] <- aic[i] + (2*n.par*(n.par+1))/(window-n.par-1)
        bic[i] <- BIC(m)
        mse[i] <- mean((m$residuals)^2)
        mm <- summary(m)
        if (all(is.finite(mm$coefficients[,1])) && all(is.finite(mm$coefficients[,4])))
          {
            coeff <- rbind(coeff,mm$coefficients[,1])
            pval <- rbind(pval,mm$coefficients[,4])
          }
        else
          {
            coeff <- rbind(coeff,rep.int(0,ncol(coeff)))
            pval <- rbind(pval,rep.int(1,ncol(pval)))
          }
      }
  }

if (! is.null(x)) 
  {
    if (c==TRUE)
      {
        colnames(coeff) <- c("const",colnames(x))
      }
    else
      {
        colnames(coeff) <- colnames(x)
      }
  }
else
  {
    coeff <- matrix(coeff,ncol=1,byrow=TRUE)
    pval <- matrix(pval,ncol=1,byrow=TRUE)
    colnames(coeff) <- c("const")
  }
colnames(pval) <- colnames(coeff)

if (is.matrix(y) || is.xts(y)) { names(y.roll.ols) <- rownames(y) } 

r <- list(y.roll.ols, aic, aicc, bic, mse, as.matrix(coeff[-1,,drop=FALSE]), as.matrix(pval[-1,,drop=FALSE]), window, as.matrix(y))
names(r) <- c("y.hat","AIC","AICc","BIC","MSE","coeff.","p.val.","window","y")
class(r) <- "reg"
return(r)

}
