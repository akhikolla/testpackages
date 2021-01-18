
stest <- function(data)
  {
    ### data - a matrix or xts object 

    if (missing(data)) { stop("please, specify data") }
    if (! is.matrix(data)) 
      { 
        data <- as.matrix(data)
        warning("data should be a matrix, the function tried to convert data to a matrix") 
      }

    x <- matrix(nrow=ncol(data),ncol=6)
    colnames(x) <- c("ADF stat.","ADF p-val.","PP stat.","PP p-val.","KPSS stat.","KPSS p-val.")
    rownames(x) <- colnames(data)
    for (i in 1:ncol(data))
      {
        x[i,1] <- adf.test(data[,i])$statistic
        x[i,2] <- adf.test(data[,i])$p.value
        x[i,3] <- pp.test(data[,i])$statistic
        x[i,4] <- pp.test(data[,i])$p.value
        x[i,5] <- kpss.test(data[,i])$statistic
        x[i,6] <- kpss.test(data[,i])$p.value
      }
    return(round(x,digits=4))
  }
  