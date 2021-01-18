
hmdmtest <- function(y,f)
  {
    ### Modified Diebold-Mariano test for one-step ahead forecasts (when presence of ARCH effects is suspected)
   
    ### a wrapper for dm.test() from "forecast" package
    
    ### y - forecasted time-series
    
    ### f - a matrix or xts object, 
    ###     first row should correspond to the selected forecast (predicted values),
    ###     next rows should correspond to the alternative forecasts (predicted values) 
    

    if (missing(f)) { stop("please, specify f") }
    if (! is.matrix(f)) 
      { 
        f <- as.matrix(f)
        warning("f should be a matrix, the function tried to convert f to a matrix") 
      }

    y <- as.vector(y)
    
    x <- matrix(nrow=(nrow(f)-1),ncol=4)
    colnames(x) <- c("HMDM stat.","HMDM p-val. different", "HMDM p-val. greater","HMDM p-val. less")    
    
    e1 <- as.vector(y) - as.vector(f[1,])
    for (i in 1:(nrow(f)-1))
      {
        e2 <- as.vector(y) - as.vector(f[i+1,])
        m <- dm.test(e1=e1,e2=e2,alternative="two.sided",h=1,power=2)$statistic
        t <- length(e1)
        m <- m / ((1 - (1/t))^0.5)
        modh <- 1 + trunc(0.5 * (t^(1/3)))
        m < - m * t^(-0.5) * (t + 1 - 2 * modh + (1/t) * modh * (modh-1))^0.5
        p.val1 <- 2 * min(pt(q=m,df=(t-1),lower.tail=FALSE), 1 - pt(q=m,df=(t-1),lower.tail=FALSE)) 
        p.val2 <- pt(q=m,df=(t-1),lower.tail=FALSE) 
        p.val3 <- 1 - pt(q=m,df=(t-1),lower.tail=FALSE) 
        
        x[i,1] <- m
        x[i,2] <- p.val1
        x[i,3] <- p.val2
        x[i,4] <- p.val3
      }
    
    return(format(round(x,digits=4),scientific=FALSE))
  }
  