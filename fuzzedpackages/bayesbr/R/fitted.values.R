#'@title Fitted Values for Theta on Beta Regression
#'@aliases fitted.values
#'@name fitted.values
#'@description A function that receives information from an estimated model uses data from the estimated theta for each iteration and returns the average of each theta in the sample.
#'@usage fitted.values(x)
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@return A vector with the average of theta estimates in the iterations (excluding warmup). The vector size is equal to the number of model observations.
#'@seealso \code{\link{bayesbr}},\code{\link{predict.bayesbr}}
fitted.values = function(x){
  n = x$info$n
  theta = x$info$samples$theta
  fitted = c()
  if(length(theta)==1){
    v_theta = theta$theta
    meanv_theta = c(round(mean(v_theta),5))
    names(meanv_theta) = "theta"
    return(meanv_theta)
  }
  else{
    for (i in 1:n) {
      v_theta = c()
      aux = paste0('theta[',i,']')
      v_theta = theta[[aux]]
      fitted = c(fitted,round(mean(v_theta),5))
    }
  }
  names(fitted) = 1:n
  return(fitted)
}
