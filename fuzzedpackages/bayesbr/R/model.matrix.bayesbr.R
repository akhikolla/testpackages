#'@title Model Matrix/Frame with All Variables for \code{bayesbr} Objects
#'@name model_matrix
#'@aliases model_matrix
#'@description The function receives all variables and their respective names, and concatenates them in a matrix.
#'@usage
#'model_matrix(object,...)
#'@param object an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param ... further arguments passed to or from other methods.
#'@return A matrix or Frame containing all variables in the model and their names used as column names.
#'@seealso \code{\link{values}},\code{\link{bayesbr}}
#'@examples
#'data("bodyfat",package="bayesbr")
#'\dontshow{
#' lines = sample(1:251,50)
#' bodyfat = bodyfat[lines,]
#' }
#'bbr = bayesbr(siri ~ wrist +I(age/100)|chest, data = bodyfat,
#'              iter = 100)
#'model_matrix(bbr)
#'
#'model_frame(bbr)
#'@export
model_matrix = function(object,...){
return(as.matrix(object$model))
}
#'@title Model Matrix/Frame with All Variables for \code{bayesbr} Objects
#'@name model_frame
#'@aliases model_frame
#'@description The function receives all variables and their respective names, and concatenates them in a matrix.
#'@usage
#'model_frame(object,...)
#'@param object an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param ... further arguments passed to or from other methods.
#'@return A matrix or Frame containing all variables in the model and their names used as column names.
#'@seealso \code{\link{values}},\code{\link{bayesbr}}
#'@examples
#'data("bodyfat",package="bayesbr")
#'\dontshow{
#' lines = sample(1:251,50)
#' bodyfat = bodyfat[lines,]
#' }
#'bbr = bayesbr(siri ~ wrist +I(age/100)|chest, data = bodyfat,
#'              iter = 100)
#'model_matrix(bbr)
#'
#'model_frame(bbr)
#'@export
model_frame = function(object,...){
  return(object$model)
}
