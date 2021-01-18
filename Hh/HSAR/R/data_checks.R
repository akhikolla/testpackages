# methods for using formula
# data checks etc

check_formula<-function(formula, data){
  
  frame <- try(suppressWarnings(model.frame(formula, data= data, na.action=na.pass)), silent=TRUE)
  if(class(frame)=="try-error") stop("The formula contains an error.", call.=FALSE)
  
  return (frame)
}

get_X_from_frame<-function(frame){
  
  X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
  if(class(X)[1]=="try-error") stop("The covariate matrix contains inappropriate values.", call.=FALSE)
  
  return(X)
}

get_y_from_frame<-function(frame){
  
  y <- model.response(frame)
  return(y)
}

check_matrix_dimensions<-function(M, n, message){
  
  if(nrow(M)!=n | ncol(M)!=n)
  {
    stop(message, call. = FALSE)
  }
}