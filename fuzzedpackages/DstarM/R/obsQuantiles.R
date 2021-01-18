#' Calculate model fit
#'
#' @param data A dataframe with: a column named \code{rt} containing response times in ms,
#' a column named \code{response} containing at most 2 response options, and an
#' optional column named \code{condition} containing a numeric index as to which conditions
#' observations belong.
#' @param probs vector of probabilities for which the corresponding values should be called
#' @param what Character. \code{'cr'} if the quantiles are to be calculated per condition-response
#' pair, \code{'c'} if the quantiles are to be calculated per condition, and
#' \code{'r'} if the quantiles are to be calculated per response.
#'
#' @description This function is nothing but a wrapper for \code{\link{quantile}}.
#'
#' @examples
#'tt = seq(0, 5, .01)
#'pars = c(.8, 2, .5, .5, .5, # condition 1
#'         .8, 3, .5, .5, .5,  # condition 2
#'         .8, 4, .5, .5, .5)  # condition 3
#'pdfND = dbeta(tt, 10, 30)
#'# simulate data
#'data = simData(n = 3e3, pars = pars, tt = tt, pdfND = pdfND)
#'probs = seq(0, 1, .01)
#'q = obsQuantiles(data, probs = probs)
#'matplot(probs, q, type = 'l', las = 1, bty = 'n')

#' @export
obsQuantiles <- function(data, probs = seq(0, 1, 0.01), what = "cr") {
  # error handling
  stopifnot(c("rt", "response") %in% names(data))
  stopifnot(is.numeric(data$rt))
  if (!(length(unique(data$response)) == 2 | length(levels(data$response)) == 
    2)) {
    stop("There need to be at least 2 response options in data$response. If only one response option has been observed, data$response should be a factor with 2 levels where the levels represent the response options.")
  }
  # calc quantiles by condition-response pair
  if (what == "cr") {
    rtList <- split(data$rt, list(data$response, data$condition))
  } else if (what == "c") {
    # by condition
    rtList <- split(data$rt, data$condition)
  } else {
    # by response
    rtList <- split(data$rt, data$response)
  }
  rtList <- lapply(rtList, stats::quantile, probs = probs)
  dd <- do.call(cbind, rtList)
  colnames(dd) <- names(rtList)
  return(dd)
}


