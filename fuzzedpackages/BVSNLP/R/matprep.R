matprep <- function(X, fixed_cols, prep, logT){
  fcols <- Filter(is.factor, X)
  allnames <- names(X)
  fnames <- names(fcols)
  numcols <- dim(fcols)[2]
  fix_flag <- as.logical(length(fixed_cols))
  fct_flag <- as.logical(numcols)
  if (fct_flag){
    outlist <- as.list(1:numcols)
    for (i in 1:numcols){
      Lv <- levels(fcols[,i])
      auxmat <- model.matrix(~fcols[,i])
      auxmat <- as.matrix(auxmat[,-1])
      colnames(auxmat) <- paste(fnames[i],Lv[-1],sep="")
      outlist[[i]] <- auxmat
    }
    fctidx <- which(allnames%in%fnames)
    if (fix_flag){
      fixfctidx <- which(fctidx%in%fixed_cols)
      if(length(fixfctidx)){
        nonfctfixidx <- fixed_cols[-which(fixed_cols%in%fctidx)]
        nf_fct_vec <- numeric(numcols)
        for (i in 1:numcols){
          Lv <- levels(fcols[,i])
          nf_fct_vec[i] <- length(Lv[-1])
        }
        mat_fct_fix <- as.matrix(as.data.frame(outlist[fixfctidx]))
        nf_fct <- sum(nf_fct_vec[fixfctidx])
        if(length(nonfctfixidx)) mat_nfct_fix <- as.matrix(X[nonfctfixidx]) else mat_nfct_fix <- NULL
        mat_nfct_nfix <- as.matrix(X[,-c(fctidx,nonfctfixidx)])
        if(length(outlist[-fixfctidx])) mat_fct_nfix <- as.matrix(as.data.frame(outlist[-fixfctidx])) else mat_fct_nfix <- NULL
        fulmat <- cbind(mat_fct_fix, mat_nfct_fix, mat_fct_nfix, mat_nfct_nfix)
        nf <- length(fixed_cols) + nf_fct - length(fixfctidx)
      } else {
        fctmat <- as.matrix(as.data.frame(outlist))
        X1 <- as.matrix(X[,fixed_cols])
        X2 <- as.matrix(X[,-c(fixed_cols,fctidx)])
        fulmat <- cbind(X1,fctmat,X2)
        nf <- length(fixed_cols)
      }
    } else {
      fulmat <- cbind(as.matrix(as.data.frame(outlist)),as.matrix(X[,-fctidx]))
      nf <- 0
    }
  } else {
    if (fix_flag){
      X1 <- X[,fixed_cols]
      X2 <- X[,-fixed_cols]
      X3 <- cbind(X1,X2);
      fulmat <- as.matrix(X3)
      nf <- length(fixed_cols)
    } else {
      fulmat <- as.matrix(X)
      nf <- 0
    }
  }
  if (prep){
    Xin <- PreProcess(fulmat, logT)
    fulmat <- Xin$X
    gnames <- Xin$gnames
  } else {
    gnames <- colnames(fulmat)
  }
  return(list(fulmat = fulmat, gnames = gnames, nf = nf))
}