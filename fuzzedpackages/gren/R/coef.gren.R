# retrieve coefficients from gren estimated model
coef.gren <- function(object, s=NULL, type=c("groupreg", "regular"), ...) {
  if(!is.null(s) & !is.numeric(s)) {
    stop("s is either NULL or a numeric")
  } else if(is.null(s)) {
    s <- object$lambda
  } else if(!any(type %in% c("groupreg", "regular"))) {
    stop("type is groupreg or regular")
  }
  
  if(type=="groupreg") {
    coefs <- coef(object$freq.model$groupreg, s=s)
  } else {
    coefs <- coef(object$freq.model$regular, s=s)
  } 
  return(as.matrix(coefs))
}