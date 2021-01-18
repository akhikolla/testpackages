#'@title Summary for \code{bayesbr} Objects
#'@aliases summary.bayesbr
#'@name summary.bayesbr
#'@description A method that receives a list of the bayesbr type and its items and displays the main information of the model, such as the residuals, a table containing statistics on the estimated coefficients and information to evaluate the quality of the model.
#'@usage \method{summary}{bayesbr}(object,type = c("","quantile","sweighted",
#'"pearson","ordinary"), prob = 0.95,...)
#'@param object an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param type A character containing the residual type returned by the model among the possibilities. The type of residue can be \emph{quantile}, \emph{sweighted}, \emph{pearson} or \emph{ordinary}. The default is \emph{quantile}.
#'@param prob a probability containing the credibility index for the HPD interval for the coefficients of the covariates.
#'@param ... further arguments passed to or from other methods.
#'@seealso \code{\link{bayesbr}},\code{\link{residuals.bayesbr}},\code{\link{print.bayesbr}},\code{\link{predict.bayesbr}}
#'@examples
#'data("bodyfat",package="bayesbr")
#'\dontshow{
#' lines = sample(1:251,15)
#' bodyfat = bodyfat[lines,]
#' }
#'bbr = bayesbr(siri ~ age + weight +
#'               wrist | biceps + forearm,
#'               data = bodyfat, iter = 100)
#'
#'summary(bbr)
#'summary(bbr, type="pearson")
#'summary(bbr, prob = 0.9)
#'summary(bbr, prob = 0.99, resid.type="sweighted")
#'\donttest{
#'bbr2 = bayesbr(siri ~ age + weight + height +
#'            wrist | biceps + forearm, data = bodyfat,
#'            iter = 100,mean_betas = 3,
#'            variance_betas = 10)
#'
#'summary(bbr2)
#'summary(bbr2, type="sweighted")
#'summary(bbr2, prob = 0.96)
#'summary(bbr2, prob = 0.95, resid.type="quantile")
#'}
#'@export
summary.bayesbr = function(object,type=c("","quantile","sweighted", "pearson","ordinary"),
                           prob=0.95,...){
  type = match.arg(type)
  if(type!=""){
    res = residuals.bayesbr(object,type)
    object$residuals = res
    object$residuals.type = type
  }
  p = object$info$p
  q = object$info$q

  cat("\nCall:", deparse(object$call, width.cutoff = floor(getOption("width") * 0.85)),
      "", sep = "\n")
  cat("Number of Chains: ",deparse(object$info$chain),"\nIter: ",deparse(object$info$iter),
      "\nWarmup:",deparse(object$info$warmup),"\n")
  types = c("quantile","sweighted", "pearson","ordinary")
  texts = c("Quantile residuals", "Standardized weighted residuals",
            "Pearson residuals","Ordinary Residuals")

  if(p>0){
  cat(sprintf("\n%s:\n", texts[types==object$residuals.type]))
  print(structure(round(as.vector(quantile(object$residuals)), digits = 5),
                  .Names = c("Min", "1Q", "Median", "3Q", "Max")))
  }
  mean = object$coefficients$mean
  if(length(mean)==0){
    cat("\nNo coefficients (in mean model)\n")
  }
  else{
    cat("\nCoefficients (mean model): \n")
    list = summary_mean(object,prob)
    table = list$table
    print.default(table)
  }
  cat("\nCoefficients (precision model): \n")
  list = summary_precision(object,prob)
  table = list$table
  print.default(table)


  if(!is.null(object$PMSE)){
    cat("\n\nPMSE: ",format(object$PMSE, digits = 5))
  }
  else{
    cat("\n")
  }
  if(p>0){
  cat("\nPseudo R-squared: ",format(object$pseudo.r.squared, digits = 5))
  cat("\nAIC: ",format(object$AIC, digits = 5))
  cat("\nBIC: ",format(object$BIC, digits = 5))
  cat("\nDIC: ",format(object$DIC, digits = 5))
  cat("\nWAIC (SE): ",format(as.vector(object$WAIC[1]), digits = 5),paste0("(",format(as.vector(object$WAIC[2]), digits = 5),")"))
  cat("\nLOOIC (SE): ",format(as.vector(object$LOOIC[1]), digits = 5),paste0("(",format(as.vector(object$LOOIC[2]), digits = 5),")"))
  }
}
