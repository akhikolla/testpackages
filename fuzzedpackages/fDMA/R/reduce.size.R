
reduce.size <- function(dma.object)
  {
    f.red <- function(x)
      {
        x$models <- x$models[1,,drop=FALSE]
        x$models[,] <- NA
        x$post.mod <- NA
        x$yhat.all.mods <- NA
        x$p.dens. <- NA
        
        return(x)
      }
    
    if (class(dma.object)=="dma")
      {
        return(f.red(dma.object))
      }
    
    if (class(dma.object)=="grid.dma")
      {
        for (i in 1:length(dma.object$models))
          {
            for (j in 1:length(dma.object$models[[i]]))
              {
                dma.object$models[[i]][[j]] <- f.red(dma.object$models[[i]][[j]])
              }
          }
        return(dma.object)
      }
  }
