
summary.altf3 <- function(object, ...)
  {
   x <- object
   cat("Mean coefficients: ")
   cat("\n")
   c <- colMeans(x$coeff.[[1]],na.rm=TRUE)
   c <- round(c,digits=4)
   print(c,quote=FALSE)
   
   cat("\n")
   cat("Frequency when p-values for t-test are less than: ")
   cat("\n")
   vv <- rep.int(NA,3)
   for (i in 1:ncol(x$p.val.[[1]]))
    {
       v <- as.vector(na.exclude(x$p.val.[[1]][,i,drop=FALSE]))
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
   colnames(vv) <- colnames(x$p.val.[[1]])
   rownames(vv) <- c("0.01","0.05","0.10")
   print(vv,quote=FALSE)
   cat("\n")

   cat("\n")
   cat("Forecast quality measures: ")
   cat("\n")
   print(x$summary,quote=FALSE)
   cat("\n")
  }
