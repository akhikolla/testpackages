#'@title Bayesian Information Criterion
#'@name BIC_bayesbr
#'@aliases BIC_bayesbr
#'@usage BIC_bayesbr(x)
#'@description A function that receives data from the estimated model, uses the information from the loglik, the number of observations of the model and the number of estimated parameters and returns the BIC, an estimator for the quality of the estimation of a model.
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@details Proposed by Stone (1979) the BIC (Bayesian Information Criterion) measures the quality of the adjustment made by the model, when comparing adjusted models with the same data, the smaller the BIC the better the adjustment.
#'
#'The BIC theory requires that the log-likelihood has been maximized, but as we are in the context of Bayesian statistics, the log-likelihood as explained in the \code{\link{logLik.bayesbr}} is made with the average of the a priori distribution for each theta and applying this value in the formula to calculate the loglik.
#'
#'The BIC is calculated by
#'\deqn{ BIC = log(n)*k - 2 * L ,}
#'
#'where \code{n} is the number of observations of the model variables, \code{k} is the number of covariates used in the model, and L is the average of the loglik chain returned by the function \code{\link{logLik.bayesbr}}.
#'@return A number corresponding to the BIC (Bayesian Information Criterion) of the estimated model.
#'@references
#'  \doi{10.1214/aos/1176344136}Schwarz, G. (1978). Estimating the dimension of a model. \emph{The annals of statistics}, \bold{6}(2), 461-464.
#'@seealso \code{\link{bayesbr}}, \code{\link{AIC_bayesbr}}, \code{\link{DIC_bayesbr}}
#'@examples
#'data("CarTask",package = "bayesbr")
#'
#'car_bayesbr <- bayesbr(probability ~ NFCCscale + task, data = CarTask,
#'                       iter =100)
#'bic = BIC_bayesbr(car_bayesbr)
#'@export
BIC_bayesbr = function(x){
  loglik = mean(x$loglik)
  p = x$info$p
  q = x$info$q
  if(p!=0){p=p-1}
  if(q!=0){q=q-1}
  k = (x$info$p) + (x$info$q)
  n = x$info$n
  BIC = log(n)*k - (2*loglik)
  return(BIC)
}
