#' Prints prediction summary of pandemic model
#'
#' S3 method designed to summarize the prediction obtained by the
#' \code{\link{posterior_predict.pandemicEstimated}} function. It is not necessary to call function
#' directly, unless it is desired to change the default arguments.
#' @method print pandemicPredicted
#' @param x An object of S3 class \code{\link{pandemicPredicted-objects}}.
#' @param summaryFun Use this function to summarize the predictions. Default argument is median, but can be any function on a vector.
#' @param printPred Valid values are 'Long' or 'Short'. Note that 'Short' will show short-term cumulative predictions, while
#' 'Long' will return the long-term daily cases predicted.
#' @param truncView How many predictions to print. Default is 3, which means that the method prints the first
#' 3 days predicted and the last 3 days of the prediction horizon.
#' @param ... Currently unused.
#' @seealso \code{\link{posterior_predict.pandemicEstimated}} and \code{\link{plot.pandemicPredicted}}.
#' @examples
#' \dontrun{
#' dataMG = load_covid("Brazil","MG")
#' estimMG = pandemic_model(dataMG)
#' predMG = posterior_predict(estimMG)
#' print(predMG, summaryFun = mean, truncView = 5)}
#' @importFrom stats median
#' @export
print.pandemicPredicted = function(x,summaryFun = stats::median,printPred = "Long",truncView = 3,...){
  cat("\nPredicted pandemic ",ifelse(x$cases_type=="confirmed", "confirmed", "death")," cases for ",x$location,". Can be plotted with plot().\n",sep="")
  if (printPred == "Long"){
    dates = max(x$data$date) + 1:ncol(x$predictive_Long)
    preds = apply(x$predictive_Long,2,summaryFun)
    names(preds) = dates
  }
  else if (printPred == "Short"){
    dates = max(x$data$date) + 1:ncol(x$predictive_Short)
    preds = apply(x$predictive_Short,2,summaryFun)
    names(preds) = dates
  }
  else stop("printPred must be \'Long\' or \'Short\'")
  cat("\nShowing predictive ",deparse(substitute(summaryFun))," for the ",tolower(printPred),"-term predictions for ",x$location,".\n\n",sep="")
  if (truncView >= (length(preds)/2) | (!truncView))
    print.default(preds)
  else {
    print.default(preds[1:truncView])
    cat("   ...\n\n")
    print.default(preds[(length(preds)-truncView+1):length(preds)])
  }
  cat("\n*For customized view, see help(print.pandemicPredicted)",sep="")
  cat("\n**For more details, see help(pandemicPredicted-xs)\n",sep="")

  invisible(x)
}
