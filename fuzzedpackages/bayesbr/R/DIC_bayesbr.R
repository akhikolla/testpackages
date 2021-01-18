#'@title Deviance Information Criterion
#'@name DIC_bayesbr
#'@aliases DIC_bayesbr
#'@description A function that receives data from the estimated model, uses the information from the loglik and returns the DIC, an estimator for the quality of the estimation of a model.
#'@usage DIC_bayesbr(x)
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@details Proposed by Spiegelhalter (2002) the DIC (Deviance Information Criterion) measures the quality of the adjustment made by the model, when comparing adjusted models with the same data, the smaller the BIC the better the adjustment.
#'
#'It is particularly useful in Bayesian model selection problems where the posterior distributions of the models have been obtained by Markov chain Monte Carlo (MCMC) simulation. DIC is an asymptotic approximation as the sample size becomes large, like AIC. It is only valid when the posterior distribution is approximately multivariate normal.
#'
#'DIC is calculate using the loglik calculated from the posterior distribution of the parameters and a calculation from the average of the posterior distribution of the parameters. To see the formula visit \href{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/1467-9868.00353}{Spiegelhalter (2002)}.
#'@return A number corresponding to the DIC (Deviance Information Criterion) of the estimated model.
#'@references
#'\doi{10.1111/1467-9868.00353} Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the royal statistical society: Series b (statistical methodology)}, \bold{64}(4), 583-639.
#'@references
#'\doi{10.1111/j.1467-9574.2005.00278.x} Van Der Linde, A. (2005). DIC in variable selection. \emph{Statistica Neerlandica}, \bold{59}(1), 45-56.
#'
#'@seealso \code{\link{bayesbr}}, \code{\link{AIC_bayesbr}}, \code{\link{BIC_bayesbr}}
#'@examples
#'data("CarTask",package="bayesbr")
#'
#'car_bayesbr <- bayesbr(probability ~ NFCCscale + task, data = CarTask,
#'                       iter =100)
#'dic = DIC_bayesbr(car_bayesbr)
#'@export
DIC_bayesbr = function(x){
  loglik = mean(x$loglik)
  d = -2*loglik
  mean_d = mean(d)
  mean_vtheta = x$fitted.values
  zeta = x$info$samples$zeta
  mean_vzeta = c()
  if(length(zeta)==1){
    v_zeta = zeta$zeta
    mean_vzeta = c(round(mean(v_zeta),5))
    names(mean_vzeta) = "zeta"
  }
  else{
    for (i in 1:x$info$n) {
      v_zeta = c()
      aux = paste0('zeta[',i,']')
      v_zeta = zeta[[aux]]
      mean_vzeta = c(mean_vzeta,round(mean(v_zeta),5))
    }
  }
  A = mean_vtheta * mean_vzeta
  B = mean_vzeta-A
  loglik2 = 0
  Y = x$y
  for(i in (1:x$info$n)){
    loglik2 = loglik2 + dbeta(Y[i],shape1 = A[i],shape2 = B[i],log = TRUE)
  }
  dbarra = -2*loglik2
  pd = mean_d - dbarra
  DIC = pd + mean_d

  return(as.numeric(DIC))
}
