predmat <- function(X){
  fcols <- Filter(is.factor, X)
  allnames <- names(X)
  fnames <- names(fcols)
  numcols <- dim(fcols)[2]
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
    fulmat <- cbind(as.matrix(as.data.frame(outlist)),as.matrix(X[,-fctidx]))
  } else {
    fulmat <- as.matrix(X)
  }
  return(fulmat)
}