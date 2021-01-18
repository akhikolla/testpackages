coef.cv.BTLLasso <- function(object, ...){
    object$coefs.repar[which.min(object$criterion),]
}

coef.BTLLasso <- function(object, ...){
  object$coefs.repar
}


logLik.cv.BTLLasso <- function(object, ...){
  ll <- object$logLik[which.min(object$criterion)]
  class(ll) <- "logLik"
  attr(ll, "df") <- object$df[which.min(object$criterion)]
  ll
}

