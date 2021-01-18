#' Plot CUSUM chart for a cusum Object
#'
#' Produces a CUSUM chart.
#' @export
#' @import utils
#' @param x An object of class cusum
#' @param signal Logical. If TRUE, signals are plotted (default)
#' @method plot cusum
#' @usage ## S3 method for class 'cusum' plot(x, signal = TRUE, ...)


plot.cusum <- function(x, signal = TRUE, ...) {
  
  assert_logical(signal, len = 1)

  if (mean(x$ct) < 0){
    y_lim <- c(min(x$ct), 0)
  } else {
    y_lim <- c(0, max(x$ct))
  }

  plot(
    x = x$t,
    y = x$ct,
    type = "n",
    xlim = c(0, max(x$t)),
    ylim = y_lim,
    ylab = expression(CUSUM[t]), xlab = "t",
    ...
  )

  if (signal == TRUE) {
    #p_signal <- which(x$signal == 1)  - 1
    p_signal <- which(x$signal == 1)
    points(
      x = x$t[p_signal],
      y = x$ct[p_signal],
      col = "orange",
      cex = 1.8,
      pch = 8,
      lwd = 2.5 
    )
  }
  
  lines(
    x = x$t,
    y = x$ct
  )

  lines(
    x = x$t,
    y = x$limit,
    col = "Blue"
  )
  
  points(
    x = x$t,
    y = x$ct,
    cex = .5,
    pch = 16
  )
}
