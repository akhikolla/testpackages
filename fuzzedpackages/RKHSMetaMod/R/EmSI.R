SI_emp<-function(res,ErrPred=NULL){
  if(!is.null(ErrPred)){
    ind <- which(ErrPred == min(ErrPred))
    var.e <- apply(res[[ind]]$`Meta-Model`$fit.v,2,var)
    SImp <- var.e/sum(var.e)
  }else if(is.null(ErrPred)){
    l <- length(res)
    SImp <- list()
    for(i in 1:l){
      var.e <- apply(res[[i]]$`Meta-Model`$fit.v,2,var)
      SImp[[i]] <- var.e/sum(var.e)
    }
  }
  return(SImp)
}#End SI_emp
