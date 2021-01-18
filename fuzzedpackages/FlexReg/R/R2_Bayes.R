#' @title Bayesian R-squared for flexreg Objects
#'
#' @description Bayesian version of R-squared for flexible regression models for proportions
#'
#' @param model an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}}.
#'
#' @details  The function provides a Bayesian version of the R-squared measure, defined as the variance of the predicted values divided by itself plus the expected variance of the errors.
#' @references {
#' Andrew Gelman, Ben Goodrich, Jonah Gabry & Aki Vehtari (2019) R-squared for Bayesian Regression Models, The American Statistician, 73:3, 307-309, DOI: 10.1080/00031305.2018.1549100
#' }
#'
#' @import rstan
#'
#' @examples
#' data("Reading")
#' FB <- flexreg(accuracy ~ iq, Reading, type="FB", n.iter=1000)
#' hist(R2_bayes(FB))
#'
#' @export
#'

R2_bayes <- function(model){
  y <- model$response
  posterior <- model$model
  y_tilde <-  rstan::extract(posterior, pars="mu", permuted=T)[[1]]
  #bayes_R2(y_tilde, y) #si potrebbe semplicemente usare la funzione di rstantools
  var_fit <- apply(y_tilde, 1, var)

  # Calcolo la distribuzione dei residui per ogni unita':
  res <- apply(rbind(y,y_tilde), 2, function(x) x[-1]-x[1])
  # Varianza dei residui condizionatamente al parametro simulato:
  var_res <- apply(res, 1, var)
  R2 <- var_fit/(var_fit + var_res)
  return(R2)
}

