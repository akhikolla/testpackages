#'@title Model Log Likelihood for \code{bayesbr} Objects
#'@name logLik.bayesbr
#'@aliases logLik.bayesbr
#'@description A function that receives the information from the estimated model, the response variable and the theta and zeta chains and returns a vector containing loglik values for each iteration excluding warmups.
#'@usage \method{logLik}{bayesbr}(object,...)
#'@param object an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param ... further arguments passed to or from other methods.
#'@details Loglik is commonly used to measure fit quality, or to assess whether an fit has converged. The loglik is calculated using maximum likelihood, but as we are in the Bayesian context we will use the mean of the posterior distribution of the parameters, so the calculation occurs from an adaptation of the original form to the loglik.
#'@return The function returns a list with \describe{\item{loglik}{A vector with the estimated model loglik chain,}\item{matrix_loglik}{A matrix with all loglik's chain.}}
#'@references
#'\doi{10.1080/0266476042000214501} Ferrari, S.L.P., and Cribari-Neto, F. (2004). Beta Regression for Modeling Rates and Proportions. \emph{Journal of Applied Statistics}, \bold{31}(7), 799–815.
#'@seealso \code{\link{bayesbr}},\code{\link{residuals.bayesbr}},\code{\link{loglikPlot}}
#'@examples
#'data("CarTask", package = "bayesbr")

#'bbr = bayesbr(probability~task + NFCCscale, iter = 100,
#'              data=CarTask, mean_betas = c(1, 0.5,1.2))
#'loglik = bbr$loglik
#'
#'loglikPlot(loglik)
#'@export
logLik.bayesbr = function(object,...){
  Y = object$y
  iter = object$info$iter
  warmup = object$info$warmup
  iterstar = iter - warmup
  n = object$info$n
  theta = object$info$samples$theta
  zeta = object$info$samples$zeta
  loglik = rep(0,iterstar)
  A = matrix(0,iterstar,n)
  B = matrix(0,iterstar,n)
  for (i in 1:n) {
    v_theta = c()
    v_zeta = c()
    if(length(theta)==1){
      v_theta = theta$theta
    }
    else{
      aux = paste0('theta[',i,']')
      v_theta = theta[[aux]]

    }
    if(length(zeta)==1){
      v_zeta = zeta$zeta
    }
    else{
      aux = paste0('zeta[',i,']')
      v_zeta = zeta[[aux]]
    }
    A[,i] = v_theta*v_zeta
    B[,i] = v_zeta-A[,i]
  }
  matrix_loglik = matrix(0,nrow=iterstar, ncol=n)
  for(i in 1:iterstar){
    for(t in 1:n){
      matrix_loglik[i,t] = dbeta(Y[t],shape1 = A[i,t],shape2 = B[i,t],log = TRUE)
      loglik[i] = loglik[i] + matrix_loglik[i,t]
    }
  }
  return(list(loglik = loglik,matrix_loglik = matrix_loglik))
}
