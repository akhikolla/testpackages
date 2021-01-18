#' Print method for \code{pandemicStats} objects
#'
#' S3 method designed for \code{pandemicStats} objects. It displays a compact summary of the
#' 95\% credible intervals for the predictions calculated by the pandemic model.
#'
#' @method print pandemicStats
#'
#' @param x An object of S3 class \code{\link{pandemicStats-objects}}.
#' @param ... Currently unused.
#'
#' @seealso \code{\link{pandemic_stats}}
#'
#' @importFrom utils head tail
#' @export

print.pandemicStats <-function(x, ...){

  cat("\n95% Credible Intervals for ",ifelse(x$data$case_type=="confirmed", "confirmed", "death"),
      " cases in ",x$data$location,"\n", sep = "")

  ST = nrow(x$ST_predict)
  LT = nrow(x$LT_predict)

  cat("\nShort-term Predictions:\n")
  width<-max(sapply(x$ST_predict, nchar))
  print(utils::head(format(x$ST_predict, width=width, justify = "right"), n=3))
  cat("   ...")
  colnames(x$ST_predict)<- NULL
  print(utils::tail(format(x$ST_predict, width=width, justify = "right"), n=3))

  cat("\n------")
  cat("\nLong-term Predictions:\n")
  width1<-max(sapply(x$LT_predict, nchar))
  print(utils::head(format(x$LT_predict, width=width1, justify = "right"), n=3))
  cat("   ...")
  colnames(x$LT_predict)<- NULL
  print(utils::tail(format(x$LT_predict, width=width1-1, justify = "right"), n=3))

  cat("\n------")
  cat("\nTotal Number of Cases:\n")
  total = data.frame(x$LT_summary$total_cases_LB, x$LT_summary$total_cases_med,
                            x$LT_summary$total_cases_UB)
  colnames(total)<- c("q2.5","med","q97.5")
  print(total, row.names = FALSE)

  cat("\n------")
  cat("\nPeak Dates:\n")
  peak = data.frame(x$LT_summary$peak_date_LB, x$LT_summary$peak_date_med,
                         x$LT_summary$peak_date_UB)
  colnames(peak)<- c("q2.5","med","q97.5")
  print(peak, row.names = FALSE)

  cat("\n------")
  cat("\nEnd Dates:\n")
  end = data.frame(x$LT_summary$end_date_LB, x$LT_summary$end_date_med,
                    x$LT_summary$end_date_UB)
  colnames(end)<- c("q2.5","med","q97.5")
  print(end, row.names = FALSE)

  cat("\n------\n")
  cat("*Use plot() to see these statistics in a graph format.\n")
  cat("*For more information, see ?'pandemicStats-xs'.\n")
  cat("*For details on the calculations, see ?pandemic_stats.\n")


  invisible(x)

}
