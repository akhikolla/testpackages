#'@title Print for \code{bayesbr} Objects
#'@name print.bayesbr
#'@aliases print.bayesbr
#'@description A method that receives a list of the bayesbr type and its items and displays the estimated coefficients.
#'@usage \method{print}{bayesbr}(x,...)
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param ... further arguments passed to or from other methods.
#'@seealso \code{\link{bayesbr}}, \code{\link{summary.bayesbr}}, \code{\link{residuals.bayesbr}}
#'@examples
#'data("bodyfat",package="bayesbr")
#'\dontshow{
#' lines = sample(1:251,50)
#' bodyfat = bodyfat[lines,]
#' }
#'bbr = bayesbr(brozek ~ wrist + density:thigh |chest, data = bodyfat,
#'              iter = 100)
#'print(bbr)
#'@export
print.bayesbr  = function(x,...){
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  if(!is.null(x$PMSE)){
    cat(paste0("\nPMSE: ",x$PMSE,"\n"))
  }
  mean = x$coefficients$mean
  precision = x$coefficients$precision
  if(length(mean)==0){
    cat("\nNo coefficients (in mean model)\n")
  }
  else{
    cat("\nCoefficients (mean model): \n")
    print.default(format(mean, digits = 5), print.gap = 2, quote = FALSE)
  }
  cat("\nCoefficients (precision model): \n")
  print.default(format(precision, digits = 5), print.gap = 2, quote = FALSE)

  cat("\nNumber of Chains: ",deparse(x$info$chain),"\nIter: ",deparse(x$info$iter),"\nWarmup:",deparse(x$info$warmup),"\n")

}
