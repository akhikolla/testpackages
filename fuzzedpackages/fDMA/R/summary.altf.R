
summary.altf <- function(object, ...)
  {
   x <- object
   
   nmods <- length(x$coeff.)
   
   if (! all(names(x$coeff.) %in% c("naive")))
    {
      cat("Mean coefficients: ")
      cat("\n")
    }
   
   c <- NULL
   n.c <- 0
   n.i <- 0
   m1 <- vector()
   
   for (i in 1:nmods)
    {
      if (names(x$coeff.)[i] %in% c("OLS","rec. OLS","roll. OLS","TVP"))
        {
          n.c <- ncol(x$coeff.[[i]])
          n.i <- n.i + 1
          m1[i] <- 1
        }
      else
        {
          m1[i] <- 0
        }
    }
   if (n.c > 0)
    {
      c <- matrix(,ncol=n.c,nrow=sum(m1))
      m1.r <- which(m1==1)
      for (i in 1:n.i)
        {
          j.r <- m1.r[i]
          c[i,] <- colMeans(x$coeff.[[j.r]],na.rm=TRUE)
        }
      j.r <- m1.r[1]
      colnames(c) <- colnames(x$coeff.[[j.r]])
      rownames(c) <- names(x$coeff.[m1.r])
        
    }
    
   if (c("MS") %in% names(x$coeff.))
    {
      j <- which(names(x$coeff.)=="MS")
      if (is.null(c))
        {
          c <- round(x$coeff.[[j]],digits=4)
          print(c,quote=FALSE)
          cat("\n")
        }
      else
        {
          c <- rbind(c,x$coeff.[[j]])
          c <- round(c,digits=4)
          print(c,quote=FALSE)
          cat("\n")
        }
    }
   else
    {
      if (! is.null(c))
        {
          c <- round(c,digits=4)
          print(c,quote=FALSE)
          cat("\n")
        }
    }


   c <- NULL
   n.c <- 0
   n.i <- 0
   m1 <- vector()
   
   for (i in 1:nmods)
    {
      if (names(x$coeff.)[i] %in% c("AR(1)","AR(2)","TVP-AR(1)","TVP-AR(2)"))
        {
          n.c <- max(ncol(x$coeff.[[i]]),n.c)
          n.i <- n.i + 1
          m1[i] <- 1
        }
      else
        {
          m1[i] <- 0
        }
    }
   if (n.c > 0)
    {
      c <- matrix(,ncol=n.c,nrow=sum(m1))
      m1.r <- which(m1==1)
      for (i in 1:n.i)
        {
          j.r <- m1.r[i]
          if (ncol(x$coeff.[[j.r]])==n.c)
            {
              c[i,] <- colMeans(x$coeff.[[j.r]],na.rm=TRUE)
            }
          else
            {
              c[i,] <- cbind(t(colMeans(x$coeff.[[j.r]],na.rm=TRUE)),NA)
            }
        }
      if (n.c==2)
        {
          colnames(c) <- c("const","ar1")
        }
      else
        {
          colnames(c) <- c("const","ar1","ar2")
        }
      rownames(c) <- names(x$coeff.[m1.r])
      c <- round(c,digits=4)
      print(c,quote=FALSE)
      cat("\n")
    }
    
    
   if (c("auto ARIMA") %in% names(x$coeff.))
    {
      j <- which(names(x$coeff.)=="auto ARIMA")
      if (is.null(c))
        {
          c <- round(x$coeff.[[j]],digits=4)
          rownames(c) <- c("auto ARIMA")
          print(c,quote=FALSE)
          cat("\n")
        }
      else
        {
          c <- round(x$coeff.[[j]],digits=4)
          rownames(c) <- c("auto ARIMA")
          print(c,quote=FALSE)
          cat("\n")
        }
    }
 


   if (! all(names(x$coeff.) %in% c("naive","TVP","TVP-AR(1)","TVP-AR(2)")))
    {
      cat("Frequency when p-values for t-test are less than: ")
    }
   cat("\n")
   
   for (k in 1:nmods)
    {
      if (! names(x$p.val.)[k] %in% c("naive","TVP","TVP-AR(1)","TVP-AR(2)","MS"))
       {
         vv <- rep.int(NA,3)
         for (i in 1:ncol(x$p.val.[[k]]))
          {
             v <- as.vector(na.exclude(x$p.val.[[k]][,i,drop=FALSE]))
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
         colnames(vv) <- colnames(x$p.val[[k]])
         rownames(vv) <- c("0.01","0.05","0.10")
         temp <- names(x$p.val)[k]
         names(temp) <- c("")
         print(temp,quote=FALSE)
         cat("\n")
         print(vv,quote=FALSE)
         cat("\n")
       }
      if (names(x$p.val.)[k] == c("MS"))
       {
         for (j in 1:2)
          {
             vv <- rep.int(NA,3)
             for (i in 1:ncol(x$p.val.[[k]]))
              {
                 v <- as.vector(na.exclude(x$p.val.[[k]][j,i,drop=FALSE]))
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
             colnames(vv) <- colnames(x$p.val[[k]])
             rownames(vv) <- c("0.01","0.05","0.10")
             temp <- paste(names(x$p.val)[k],"Regime",as.character(j))
             names(temp) <- c("")
             print(temp,quote=FALSE)
             cat("\n")
             print(vv,quote=FALSE)
             cat("\n")
          }
       }
    }
 
 
 
   cat("\n")
   cat("Forecast quality measures: ")
   cat("\n")
   print(x$summary,quote=FALSE)
   cat("\n")
  }
