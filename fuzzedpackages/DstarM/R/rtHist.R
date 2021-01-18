#' Make histograms of reaction time data
#'
#' @param data A reaction time dataset. Must be a dataframe with $rt, $condition
#' and $response.
#' @param what @param what What to plot. Can be 'cr' for 'condition-response pairs,
#' 'c' for condition, and 'r' for response.
#' @param layout An optional layout.
#' @param nms An optional vector of names for each plot. If omitted the names
#' will be based on the contents of \code{data$condition} and/or \code{data$response}.
#' @param ggplot ggplot Logical, should ggplot2 be used instead of base R graphics? If set to TRUE,
#' some arguments from \code{linesArgs} and \code{...} will be ignored (but can be added
#' to plots manually).
#' @param ... Arguments to be passed to \code{hist}
#'
#' @return invisible()
#'
#' @details This function and \code{\link{rtDescriptives}} are helper functions to inspect raw data.
#'
#' @examples
#' tt = seq(0, 5, .01)
#' dat = simData(n = 3e4, pars = rep(.5, 5), tt = tt, pdfND = dbeta(tt, 10, 30))
#' rtHist(dat, breaks = tt, xlim = c(0, 1))

#' @export
rtHist <- function(data, what = "cr", layout = NULL, nms = NULL, ggplot = FALSE, 
  ...) {
  rt <- switch(what, cr = split(data$rt, list(data$response, data$condition)), 
    c = split(data$rt, data$condition), r = split(data$rt, data$response))
  if (is.null(layout)) 
    layout <- matrix(1:length(rt), nrow = floor(sqrt(length(rt)))) else stopifnot(length(unique(c(layout))) == length(rt))
  if (is.null(nms)) 
    nms <- names(rt) else stopifnot(length(nms) == length(rt))
  
  histArgs <- list(...)
  histDefArgs <- list(xlab = "Reaction Time", las = 1, bty = "n")
  idx <- unlist(lapply(histArgs[names(histDefArgs)], is.null), use.names = FALSE)
  histArgs[names(histDefArgs)[idx]] <- histDefArgs[idx]
  
  if (!ggplot) {
    layout(layout)
    for (i in 1:length(rt)) {
      histArgs$x <- rt[[i]]
      histArgs$main <- nms[i]
      do.call(graphics::hist, histArgs)
    }
    layout(1)
    return(invisible())
  } else {
    plotList <- vector("list", length(rt))
    for (i in 1:length(rt)) {
      histDat <- data.frame(rt = rt[[i]])
      plotList[[i]] <- ggplot2::ggplot(data = histDat, ggplot2::aes_string(x = "rt")) + 
        ggplot2::geom_histogram(ggplot2::aes_string(y = "..density.."), 
          breaks = histArgs$breaks) + ggplot2::scale_x_continuous(name = histArgs$xlab, 
        limits = histArgs$xlim) + ggplot2::scale_y_continuous(name = histArgs$ylab) + 
        ggplot2::ggtitle(nms[i])
    }
    return(plotList)
  }
}


