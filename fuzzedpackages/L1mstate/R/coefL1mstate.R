coefl1mstate=function(object,s=c("lambda.pcvl","lambda.min")){
  if(is.numeric(s)){
    check = length(which(object$fit[,1]==s))
    if(check!=0){
      ind=which(object$fit[,1]==s)
      return(object$aBetaO[[ind]])
    }else{
      stop("Invalid value of s")
    }
  }else{
    if(is.character(s)){
      s=match.arg(s)
      if(s=="lambda.pcvl") return(object$pBetaO)
      if(s=="lambda.min") return(object$mBetaO)
    }else stop("Invalid form for s")
  } 
}
