
archtest <- function(ts,lag=NULL)
  {
    ### Lagrange Multiplier (LM) test for the presence of ARCH effects
   
    ### ts - vector, tested time-series
    ### lag - suspected order of ARCH process
    
    if (missing(ts)) { stop("please, specify ts") }
    
    n <- deparse(substitute(ts))
    
    if (! is.vector(ts)) 
      { 
        ts <- as.vector(ts)
        warning("ts should be a vector, the function tried to convert ts to a vector") 
      }
    if (is.null(lag)) { lag <- 1 }

    e <- (lm(ts ~ 1)$residuals)^2
    e.reg <- as.vector(e)
    
    f.lags <- function(v,l)
      {
        v.na <- rep(NA,l)
        lag.v <- c(v.na,v)
        lag.v <- lag.v[1:(length(lag.v)-l)]
        return(lag.v)
      }
    
    for (i in 1:lag)
      {
        e.reg <- rbind(e.reg,f.lags(e,i))
      }
    
    e.reg <- e.reg[,-c(1:lag)]
    e.reg <- t(e.reg)
    
    m <- summary(lm(e.reg[,1] ~ e.reg[,-1]))
    st <- m$r.squared * (length(ts)-lag)
    p.val <- pchisq(q=st,df=lag,lower.tail=FALSE)
    
    names(st) <- "statistic"
    names(lag) <- "lag"
    ret <- list(st,lag,paste("ARCH effects of order",lag,"are present"),p.val,"Engle's LM ARCH Test",n)
    names(ret) <- c("statistic","parameter","alternative","p.value","method","data.name")
    class(ret) <-"htest"
    
    return(ret)
  }
  