
summary.dma <- function(object, ...)
  {
   
   x <- object
   
   E.coef <- round(colMeans(x$exp.coef.),digits=4)
   probs <- round(colMeans(x$post.incl),digits=2)
   min.probs <- round(apply(x$post.incl,2,min),digits=2)
   max.probs <- round(apply(x$post.incl,2,max),digits=2)
   
   inc <- vector()
   
   for (i in 1:length(probs))
    {
      inc[i] <- round(length(x$post.incl[x$post.incl[,i]>0.5,i]) / nrow(x$post.incl),digits=2)
    }
   rm(i)
   
   s <- cbind(E.coef,min.probs,probs,max.probs,inc)
   colnames(s) <- c("E(coeff.)","min pip", "mean pip", "max pip", "inc")
   rownames(s) <- colnames(x$models)
   if (! x$parameters[1,4] == "DMA") { s <- s[,c(1,5)] }
   
   cat("Model: ")
   cat(colnames(x$y))
   cat(" ~ ")
   cat(paste(colnames(x$models),collapse=" + "))
   cat("\n")
   cat("\n")
     
   print(x$parameters,quote=FALSE)
   cat("\n")
   
   err <- rbind(x$RMSE,x$MAE)
   colnames(err) <- c("model")
   err <- cbind(err,x$benchmarks)
   print(round(err,digits=4))

   cat("\n")
   cat("observations: ")
   cat(length(x$y.hat))
   cat("\n")
   cat("models: ")
   cat(nrow(x$models))
   cat("\n")
   cat("variables (incl. constant): ")
   cat(ncol(x$models))
   
   cat("\n")
   cat("\n")
   print(s,quote=FALSE)
   cat("\n")
   
   if (x$parameters[1,4] == "DMA") 
    { 
      cat("* pip = posterior inclusion prob.")
      cat("\n")
      cat("* inc = relative frequency of a variable posterior inclusion prob. > 1/2")
    }
   if (x$parameters[1,4] == "DMS" || x$parameters[1,4] == "MED")
    {
      cat("* inc = relative frequency of a variable inclusion")
    }
   cat("\n")
  }
