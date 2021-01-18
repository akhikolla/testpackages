#' @title keras_r2
#'
#' @description This function is to calculate rsquared value for regression models.
#'
#' @param y_true tensor of ground truth   
#'
#' @param  y_pred tensor of predicted values
#'
#' @return keras tensor, double for rsquared
#'
#' @export keras_r2
keras_r2=function(y_true, y_pred) {
  SS_res =keras::k_sum(keras::k_square(y_true-y_pred ))
  SS_tot =keras::k_sum(keras::k_square(( y_true - keras::k_mean(y_true))))
  return ( 1 - SS_res/(SS_tot + keras::k_epsilon()))
}

#' @title r2_squ
#'
#' @description This function is to calculate rsquared value for regression models.
#'
#' @param obs vector, observed values   
#'
#' @param  res residual (observed values-predicted values)
#'
#' @return double for rsquared
#'
#' @export r2_squ
r2_squ=function (obs, res){
  yy = obs - matrix(mean(obs), nrow = nrow(array(obs)))
  r2 = 1 - (t(res) %*% res)/(t(yy) %*% yy)
  return(r2[1,1])
}

#' @title rmse
#'
#' @description This function is to calculate rmse value for regression models.
#'
#' @param obs vector, observed values   
#'
#' @param  pre vector, predicted values
#'
#' @return double for rmse
#'
#' @export rmse
rmse=function (obs, pre) {
  error = obs - pre
  ret = sqrt(mean(error^2))
  return(ret)
}

#' @title rSquared
#'
#' @description This function is to calculate rsquared value for regression models.
#'
#' @param obs vector, observed values   
#'
#' @param  res residual (observed values-predicted values)
#'
#' @return double for rSquared
#'
#' @export r2
rSquared=function (obs, res){
  yy = obs - matrix(mean(obs), nrow = nrow(array(obs)))
  r2 = 1 - (t(res) %*% res)/(t(yy) %*% yy)
  return(r2)
}

