# function to make predictions with the gren object
predict.gren <- function(object, newx, unpenalized=NULL, s=NULL, 
                         type=c("groupreg", "regular"), ...) {
  if(class(object)!="gren") {
    stop("object should be a gren fit")
  } else if(!is.null(s) & !is.numeric(s)) {
    stop("s is either NULL or a numeric")
  } else if(is.null(s)) {
    s <- object$lambda
  } else if(!any(type %in% c("groupreg", "regular"))) {
    stop("type is either groupreg or regular")
  }
  
  if(is.data.frame(unpenalized)) {
    unpenalized <- model.matrix(as.formula(paste("~", paste(colnames(
      unpenalized), collapse="+"))), data=unpenalized)[, -1]
  }
  x <- cbind(unpenalized, newx)
  if(type=="groupreg") {
    prob <- predict(object$freq.model$groupreg, x, s=s, type="response")
  } else {
    prob <- predict(object$freq.model$regular, x, s=s, type="response")
  }
  return(prob)
}

