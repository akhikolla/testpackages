mdr <- function(X, status, fold=2, t=NULL, top=3, NAvalues=NA, cvc=0, fix=NULL, verbose=FALSE){
  
    if(is.data.frame(X)){
      X <- as.matrix(X)
      if(verbose) message("X casted from data.frame to matrix.")
    } 
    if(!is.matrix(X)) stop("Matrix required as input for X")
    
    if(is.na(NAvalues)){
      X[is.na(X)] <- 9999    # Use an unlikly placeholder for missing data
    } else {
      if(sum(is.na(X))>0 ) stop("NA values in snpmatrix. Please use NAvalues=NA as option.")
      X[X==NAvalues] <- 9999
    }
   if(verbose==TRUE) message("Missing data detected: ", sum(X==NAvalues))
  
    N <- nrow(X)
    
    if(is.character(fix)){
      fix <- which((colnames(X)==fix)==TRUE)
    }

    ifelse(is.null(fix), fix <- -1, fix <- fix - 1)
    
    if(is.null(t)) t <- table(status)[2]/table(status)[1]
    if(nrow(X)!=length(status)) stop ("nrow(X) / length(status) mismatch!\n")

    res <- mdr.C(X=X, fold=fold, status=status, t=t, cv=0, cvp=1, top=top, na=as.numeric(NAvalues), fix)  

    res.sampled <- list()
    if(cvc>0){
      cvcBorders <- floor(seq(0,nrow(X),length.out = cvc +1))
      for(i in 1:cvc){
        X.sampled <- X[(cvcBorders[i]+1):cvcBorders[i+1],]
        status.sampled <- status[(cvcBorders[i]+1):cvcBorders[i+1]]
        res.sampled[[i]] <- mdr.C(X=X.sampled, fold=fold, status=status.sampled, t=t, cv=0, cvp=1, top=top, na=as.numeric(NAvalues), fix)  
      }
    }

    result <- list(mdr=res,fold=fold,t=t,top=top,fix=fix,X=X,status=status, res.sampled=res.sampled)

    class(result) <- "mdr"
    result
} 

