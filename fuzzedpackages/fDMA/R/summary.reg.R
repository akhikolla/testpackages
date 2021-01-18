
summary.reg <- function(object, ...)
  {
   x <- object
   cat("Mean coefficients: ")
   c <- colMeans(x$coeff.,na.rm=TRUE)
   cat("\n")
   print(round(c,digits=4),quote=FALSE)
   cat("\n")
   cat("Frequency when p-values for t-test are less than: ")
   cat("\n")
   vv <- rep.int(NA,3)
   for (i in 1:ncol(x$p.val.))
    {
       v <- as.vector(na.exclude(x$p.val.[,i,drop=FALSE]))
       v1 <- v[v<0.01]
       v2 <- v[v<0.05]
       v3 <- v[v<0.10]
       v1 <- length(v1) / length(v)
       v2 <- length(v2) / length(v)
       v3 <- length(v3) / length(v)
       v <- c(v1,v2,v3)
       vv <- rbind(vv,v)
    }
   vv <- vv[-1,,drop=FALSE]
   vv <- t(vv)
   vv <- round(vv,digits=2)
   colnames(vv) <- colnames(x$p.val.)
   rownames(vv) <- c("0.01","0.05","0.10")
   print(vv,quote=FALSE)
   cat("\n")
   e <- (mean((as.vector(x$y)-as.vector(x$y.hat))^2,na.rm=TRUE))^0.5
   e <- round(e,digits=4)
   cat("RMSE: ",e)
   cat("\n")
   e <- mean(abs(as.vector(x$y)-as.vector(x$y.hat)),na.rm=TRUE)
   e <- round(e,digits=4)
   cat("MAE:  ",e)
   cat("\n")
   if (! is.null(x$window)) { cat("rolling window: ",x$window) } 
   cat("\n")
  }
