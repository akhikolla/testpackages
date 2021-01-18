
dmtest <- function(y,f)
  {
    ### Diebold-Mariano test for one-step ahead forecasts
   
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
    colnames(x) <- c("DM stat.","DM p-val. different", "DM p-val. greater","DM p-val. less")
    
    e1 <- as.vector(y) - as.vector(f[1,])
    for (i in 1:(nrow(f)-1))
      {
        e2 <- as.vector(y) - as.vector(f[i+1,])
        m <- dm.test(e1=e1,e2=e2,alternative="two.sided",h=1,power=2)
        st <- m$statistic / ((1 - (1/length(e1)))^(1/2))
        
        x[i,1] <- st
        
        x[i,2] <- 2 * min(pnorm(st,lower.tail=FALSE),1 - pnorm(st,lower.tail=FALSE))
        x[i,3] <- pnorm(st,lower.tail=FALSE)
        x[i,4] <- 1 - pnorm(st,lower.tail=FALSE)
      }
      
    return(format(round(x,digits=4),scientific=FALSE))
  }
  