

print.qbnmix <- function(x, ...)
  {
   cat("Estimated means: ")
   c <- x$mu[[1]][nrow(x$mu[[1]]),]
   for (i in 1:length(x$mu))
     { 
      c <- rbind(c,x$mu[[i]][nrow(x$mu[[i]]),])
     }
   c <- c[-1,]
   rownames(c) <- rep(" ",nrow(c))
   colnames(c) <- rep(" ",ncol(c))
   cat("\n")
   print(format(round(c,digits=4),scientific=FALSE),quote=FALSE)
   cat("\n")
   cat("Estimated covariance matrices: ")
   cat("\n")
   for (i in 1:length(x$R))
     { 
      p <- x$R[[i]]
      colnames(p) <- rep(" ",ncol(p))
      rownames(p) <- rep(" ",nrow(p))
      print(format(round(p,digits=4),scientific=FALSE),quote=FALSE)
      cat("\n")
     }
  }

