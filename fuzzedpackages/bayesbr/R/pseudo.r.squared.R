#'@title Pseudo R Squared Calculate
#'@name pseudo.r.squared
#'@aliases pseudo.r.squared
#'@description The function receives the model information, as well as the variable response and the predicted theta values and calculates the model's pseudo.r.squared, using the formula proposed by Cribarri-Neto and Ferrari.
#'@usage pseudo.r.squared(x)
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@details Ferarri and Cribari-Neto (2004) defined the pseudo.r.squared as the square of the correlation between the theta estimated by the maximum likelihood and the logis of the variable response of the model. But as we are in the context of Bayesian statistics, the estimated theta is given by the mean of the posterior distribution of the parameter. So the informed pseudo.r.squared is a Bayesian adaptation to what was suggested by Ferarri and Cribari-Neto (2004).
#'@return A number containing the pseudo r squared of the adjusted model, this value can be used to assess the quality of the model.
#'@references
#'\doi{10.1080/0266476042000214501} Ferrari, S., & Cribari-Neto, F. (2004). Beta regression for modelling rates and proportions. \emph{Journal of applied statistics}, \bold{31}(7), 799-815.
#'@seealso \code{\link{bayesbr}},\code{\link{fitted.values}},\code{\link{AIC_bayesbr}}
#'

pseudo.r.squared = function(x){
  Y = x$y
  n = x$info$n
  theta = x$info$samples$theta
  vet_means_theta = c()
  for (i in 1:n) {
    v_theta = c()
    if(length(theta) == 1){
      v_theta = theta$theta
    }
    else{
      aux = paste0('theta[',i,']')
      v_theta = theta[[aux]]
    }
    mean_theta = mean(v_theta)
    vet_means_theta = c(vet_means_theta,mean_theta)
  }
  return(cor(qlogis(Y),vet_means_theta)^2)
}
