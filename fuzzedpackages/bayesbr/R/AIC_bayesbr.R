#'@title Akaike Information Criterion
#'@name AIC_bayesbr
#'@aliases AIC_bayesbr
#'@description A function that receives the estimated model data, uses the information from the loglik and the number of estimated parameters and returns the AIC, an estimator for the quality of the estimation of a model.
#'@usage AIC_bayesbr(x)
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@details
#'Proposed by Akaike (1974) the AIC (Akaike Information Criterion) measures the quality of the adjustment made by the model, when comparing adjusted models with the same data, the smaller the AIC the better the adjustment.
#'
#'The AIC theory requires that the log-likelihood has been maximized, but as we are in the context of Bayesian statistics, the log-likelihood as explained in the \code{\link{logLik.bayesbr}} is made with the average of the a priori distribution for each theta and applying this value in the formula to calculate the loglik.
#'The AIC is calculated by
#' \deqn{AIC = 2 * k - 2 * L ,}
#' where \code{k} is the number of covariates used in the model, and \code{L} is the average of the loglik chain returned by the function \code{\link{logLik.bayesbr}}.
#'
#' @return A number corresponding to the AIC (Akaike Information Criterion) of the estimated model.
#' @references
#' \doi{10.1109/TAC.1974.1100705} Akaike, H. (1974). A new look at the statistical model identification. \emph{IEEE transactions on automatic control}, \bold{19}(6), 716-723.
#' @seealso \code{\link{logLik.bayesbr}},\code{\link{BIC_bayesbr}},\code{\link{DIC_bayesbr}}
#' @examples
#' data("CarTask",package = "bayesbr")
#'
#' car_bayesbr <- bayesbr(probability ~ NFCCscale + task,
#'                       data = CarTask,iter =100)
#' aic = AIC_bayesbr(car_bayesbr)
#'@export
AIC_bayesbr = function(x){
  loglik = mean(x$loglik)
  p = x$info$p
  q = x$info$q
  n = x$info$n

  if(p!=0){p=p-1}
  if(q!=0){q=q-1}
  k = (x$info$p) + (x$info$q)
  AIC = log(k) - (2*loglik)
  return(AIC)
}
