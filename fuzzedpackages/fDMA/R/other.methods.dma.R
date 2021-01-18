
predict.dma <- function(object,newdata,type,...)
  {
    if (class(object)=="dma")
      {
        co <- coef.dma(object)
        newdata <- cbind(1,newdata)
        
        if (type=="backward")
          {
            pred <- as.vector(rowSums(co*newdata))
          }
        if (type=="forward")
          {
            co <- co[nrow(co),]
            pred <- vector()
            for (i in 1:nrow(newdata))
              {
                pred[i] <- sum(co*newdata[i,])
              }
          }
        
        return(pred)
      }
  }

fitted.dma <- function(object,...)
  {
    if (class(object)=="dma")
      {
        return(as.vector(object$y.hat))
      }
  }

residuals.dma <- function(object,...)
  {
    if (class(object)=="dma")
      {
        return(as.vector(object$y.hat)-as.vector(object$y))
      }
  }

coef.dma <- function(object,...)
  {
    if (class(object)=="dma")
      {
        return(object$exp.coef.)
      }
  }

rvi <- function(dma.object)
  {
    if (class(dma.object)=="dma")
      {
        return(dma.object$post.incl)
      }
  }
