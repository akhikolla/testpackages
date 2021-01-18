mModeTSIRMatrix <-
  function(x, m, y, h, ...){
    xm <- mFlatten(x, m)
    
    dim_xm <- dim(xm)
    p1 <- dim_xm[1]
    p2 <- dim_xm[2]
    n <- dim_xm[3]
    
    if(is.factor(y)){
      slices <- sapply(y, match, unique(unlist(y)))
      h <- length(levels(y))
    }
    else{
      slices <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 1, by = 1/h), ...)),
                            include.lowest = TRUE ,labels = FALSE))
    }
    
    slicemean <- sapply(1:h, function(i) apply(xm[, , which(slices == i)], 1:2, mean))
    dim(slicemean) <- c(p1, p2, h)
    
    nh <- as.numeric(table(slices))
    return(mModeCovariance(sweep(slicemean, 3, sqrt(nh/n), "*"), 1, center = TRUE)/p2)
    
  }

