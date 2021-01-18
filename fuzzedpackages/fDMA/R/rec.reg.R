
rec.reg <- function(y,x=NULL,c=NULL)
  {

    ### estimates recursive regression model

    ### y - a numeric or a column matrix of a dependent variable

    ### x - a matrix of independent variables (drivers), different columns correspond to different variables

    ### c - a parameter indicating whether constant is included,
    ###     by default c=TRUE (constant is included),
    ###     it is not possible to set c=FALSE if x=NULL 


    if (is.null(c)) { c <- TRUE }
    
    if (is.null(x)) { c <- TRUE }
    
    if (! is.null(x)) { x <- as.matrix(x) }
    y.old <- y
    y <- as.vector(y)

    fv <- vector()
    aic <- vector()
    aicc <- vector() 
    bic <- vector() 
    mse <- vector() 
    if (! is.null(x))
      {
        if (c==TRUE)
          {
            coeff <- matrix(NA,ncol=(ncol(x)+1),nrow=1)
            n.par <- ncol(x)+2
          }
        else
          {
            coeff <- matrix(NA,ncol=ncol(x),nrow=1)
            n.par <- ncol(x)+1
          }
      }
    else
      {
        coeff <- matrix(NA,ncol=1,nrow=1)
        n.par <- 2
      }
    pval <- coeff



    if (is.null(x))
      {
        m <- lm(y[1] ~ 1)
      }
    else
      {
        if (c==TRUE)
          {
            m <- lm(y[1] ~ t(x[1,]))
          }
        else
          {
            m <- lm(y[1] ~ t(x[1,]) -1)
          }
      }
    fv[1] <- m$fitted.values[1]
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

    for (i in 2:length(as.vector(y)))
      {
        if (is.null(x))
          {
            m <- lm(y[1:i] ~ 1)
            n.par <- 2
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
          }
        fv[i] <- m$fitted.values[i]
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

      if (is.matrix(y.old) || is.xts(y.old)) { names(fv) <- rownames(y.old) } 

      r <- list(fv, aic, aicc, bic, mse, as.matrix(coeff[-1,,drop=FALSE]), as.matrix(pval[-1,,drop=FALSE]), as.matrix(y.old))
      names(r) <- c("y.hat","AIC","AICc","BIC","MSE","coeff.","p.val.","y")
      class(r) <- "reg"
      return(r)

  }
