#' Check that the basic input are reasonable
#' 
#' @param fitLME : fitted Linear Mixed Effects model
#' @param fitCOX : fitted Proportional Hazards Regression model 
#' @param rho :  number specifying the logarthmic transformation model 
#' @return VOID

CheckInputs <- function(fitLME, fitCOX, rho){

  if(!inherits(fitLME, "lme"))
    stop("\n'fitLME'must be a fit returned by lme().")
  if(length(fitLME$group) > 1)
    stop("\n nested random-effects are not allowed in lme().")
  if(!is.null(fitLME$modelStruct$corStruct))
    warning("\n correlation structure in 'fitLME' is ignored.")
  if(!is.null(fitLME$modelStruct$varStruct))
    warning("\n heteroscedasticity structure in 'fitLME' is ignored.")
  
  if(!inherits(fitCOX, "coxph"))
    stop("\n'fitCOX' must be a fit returned by coxph().")
  if(is.null(fitCOX$x))
    stop("\n must specify argument 'x=TRUE' when using coxph().")
  
  if(rho < 0) {
    rho <- 0 # fit Cox model if not specified #
    warning("\n rho<0 is not valid, Cox model is fitted instead!")
  }
}
