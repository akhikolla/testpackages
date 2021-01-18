

compareModels.iModel <- function(object, object2, k=2, ...){
  comp <- compareGC(.glist(object), .glist(object2))
  logL1 <- logLik(object)
  logL2 <- logLik(object2)
  
  if (!comp["comparable"]){
    cat("Models are not nested and can hence not be compared\n")
    return()
  } else {
    smaller <- which(comp[1:2])
    if (smaller==1){
      larger <- 2
      d.logL <- c(logL2-logL1)
      d.df   <- c(attr(logL1,"df")-attr(logL2, "df"))
      large  <- .glist(object2)
      small  <- .glist(object)
    } else {
      larger <- 1
      d.logL <- c(logL1-logL2)
      d.df   <- c(attr(logL2,"df")-attr(logL1, "df"))
      large  <- .glist(object)
      small  <- .glist(object2)
    }  
  }
      
  ans <- list(small=small, large=large, d.logL=d.logL, d.df=d.df, 
              p.value=1-pchisq(d.logL, df=d.df),
              k=k, 
              aic=2*d.logL-k*d.df)      
    
  class(ans) <- "compareiModels"
  ans
}

print.compareiModels <- function(x,...){
  cat("Large:\n")
  str(x$large, give.head=FALSE,no.list=TRUE,comp.str=" ")
  cat("Small:\n")
  str(x$small, give.head=FALSE,no.list=TRUE,comp.str=" ")
  cat(sprintf("-2logL: %8.2f df: %i AIC(k=%4.1f): %8.2f p.value: %f\n",
              2*x$d.logL, x$d.df, x$k, x$aic, x$p.value))
}


compareGC <- function(gc1, gc2){

  ## Are all elements of gc1 contained in an element of gc2 ??
    ##gc.1in2 <- all(unlistPrim(lapply(gc1, function(zzz) isin(gc2, zzz))))
    gc.1in2 <- all(unlistPrim(lapply(gc1, function(zzz) is_inset(zzz, gc2))))
  ## Are all elements of gc2 contained in an element of gc1 ??
    ##gc.2in1 <- all(unlistPrim(lapply(gc2, function(zzz) isin(gc1, zzz))))
    gc.2in1 <- all(unlistPrim(lapply(gc2, function(zzz) is_inset(zzz, gc1))))

  same <- gc.1in2 & gc.2in1
  comp <- sum(c(gc.1in2, gc.2in1))>0

  ans <- c(gc.1in2=gc.1in2, gc.2in1=gc.2in1, identical=same, comparable=comp)
  ans
}

compareGCpairs <- function(gc1, gc2){
  
  gc1.pairs <- do.call(cbind, lapply(gc1, combn_prim, 2))
  gc2.pairs <- do.call(cbind, lapply(gc2, combn_prim, 2))
  
  gc1.pairs <- split(gc1.pairs, col(gc1.pairs))
  gc2.pairs <- split(gc2.pairs, col(gc2.pairs))
  
  compareGC(gc1.pairs, gc2.pairs)
}


