#'@title Matrix with All Variables for \code{bayesbr} Objects
#'@name model.bayesbr
#'@aliases model.bayesbr
#'@description The function receives all variables and their respective names, and concatenates them in a matrix.
#'@usage
#'\method{model}{bayesbr}(Y,X = NULL,W = NULL,name_y,names_x = NULL,
#'              names_w = NULL)
#'@param Y A vector containing the model response variable,
#'@param X A matrix containing the covariates for theta of the model,
#'@param W A matrix containing the covariates of the model for zeta,
#'@param name_y The name passed in the call to the \code{\link{bayesbr}} function for the variable response,
#'@param names_x The name passed in the call to the \code{\link{bayesbr}} function for the covariates for theta,
#'@param names_w The name passed in the call to the \code{\link{bayesbr}} function for covariates for zeta.
#'@return A matrix containing all variables in the model and their names used as column names.
#'@seealso \code{\link{values}},\code{\link{bayesbr}}
model.bayesbr = function(Y,X = NULL,W = NULL,name_y,names_x = NULL,names_w = NULL){
  model = Y
  verification = name_y
  if(!is.null(X)){
    if(is.null(ncol(X))){
      tam = 1
    }
    else{
      tam = ncol(X)
    }
    for (i in 1:tam) {
      aux = TRUE
      for(name in verification){
        if (name == names_x[i]) {
          aux = FALSE
          break
        }
      }
      if(aux && names_x[i]!="(Intercept)"){
        if(!is.null(ncol(X))){
          model = cbind(model,X[,i])
        }
        else{
          model = cbind(model,X[i])
        }
        verification = c(verification,names_x[i])
      }
    }
  }

  if(!is.null(W)){
    if(is.null(ncol(W))){
      tam = 1
    }
    else{
      tam = ncol(W)
    }
    for (i in 1:tam) {
      aux = TRUE
      for(name in verification){
        if (name == names_w[i]) {
          aux = FALSE
          break
        }
      }
      if(aux && names_w[i]!="(Intercept)"){
        if(!is.null(ncol(W))){
          model = cbind(model,W[,i])
        }
        else{
          model = cbind(model,W[i])
        }
        verification = c(verification,names_w[i])
      }
    }
  }
  model = data.frame(model)
  colnames(model) = verification
  return(model)
}
