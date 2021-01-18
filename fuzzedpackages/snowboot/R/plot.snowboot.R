#' Plot Degree Distribution Estimates
#'
#' Plot LSMI-based point estimates of probabilities of node degrees, \eqn{\hat{f}(k)}, and
#' of mean degree, \eqn{\hat{\mu}}, where \eqn{k = 0, 1, \ldots} are the degrees.
#' The point estimates are supplemented with box-and-whisker plots
#' of bootstrapped values (if the input is a \code{\link{boot_dd}} output) or element-wise
#' bootstrap confidence intervals (if the input is a \code{\link{boot_ci}} output).
#' See \insertCite{chen_etal_2018_snowboot;textual}{snowboot}.
#'
#' @param x output of \code{\link{lsmi_dd}}, \code{\link{boot_dd}}, or \code{\link{boot_ci}}.
#' @param k an integer vector with degrees to plot.
#' By default, all degrees represented in \code{x} are plotted.
#' @param plotmu logical value indicating whether to plot the results for mean degree
#' (default is \code{TRUE}).
#' @param plotlegend logical value indicating whether to plot a legend
#' (default is \code{TRUE}).
#' @param col0 color to plot horizontal zero-line at \eqn{f(k) = 0}.
#' Use \code{NA} for no plotting.
#' @param lwd0 width of the horizontal zero-line at \eqn{f(k) = 0}.
#' @param colpt color for plotting point estimates.
#' @param lwdpt line width for plotting point estimates.
#' @param pchpt point type for plotting point estimates
#' (see argument \code{pch} in \code{\link[graphics]{points}}).
#' @param coli color for plotting lines or borders of box-plots for bootstrap estimates.
#' @param colibg background color, if plotting boxplots of bootstrapped estimates
#' (see argument \code{border} in \code{\link[graphics]{boxplot}}).
#' @param length length of arrows, if plotting bootstrap confidence intervals
#' (see argument \code{length} in \code{\link[graphics]{arrows}}).
#' @param boxwex argument of \code{\link[graphics]{boxplot}} function.
#' @param legendargs additional arguments for plotting the legend
#' (see arguments in \code{\link[graphics]{legend}}).
#' @param las argument of \code{\link[graphics]{plot}} function.
#' @param ... additional arguments to pass to the \code{\link[graphics]{plot}} function.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#'
#' @examples
#' net <- artificial_networks[[1]]
#' x <- lsmi_dd(net = net, n.wave = 2, n.seed = 40)
#' plot(x)
#'
#' x2 <- boot_dd(x)
#' plot(x2, k = c(1:10))
#'
#' x3 <- boot_ci(x2, prob = 0.99)
#' plot(x3, k = c(1:10))
#'
plot.snowboot <- function(x, k = NULL, plotmu = TRUE,
                          plotlegend = TRUE,
                          col0 = "gray50", lwd0 = 1,
                          colpt = "royalblue3", lwdpt = 2, pchpt = 4,
                          coli = "palegreen3", colibg = "palegreen",
                          length = 0.1,
                          boxwex = 0.4,
                          legendargs = list(x = "topright", cex = 0.9, bty = "n"),
                          las = 1,
                          ...){
  if (is.null(k)) {
    k <- 0:(length(x$fk) - 1)
  }
  k <- sort(k)
  k <- k[k <= (length(x$fk) - 1)]
  if (length(k) == 0) {stop(paste("Specify k within the range from 0 to", (length(x$fk) - 1)))}
  ind <- is.element(0:(length(x$fk) - 1), k)
  YLIM <- c(0, max(x$fk[ind]))
  if (plotmu) YLIM[1] = -YLIM[2]/4
  plotbp <- plotci <- FALSE
  inputlist <- names(x)
  legtext <- c("Point estimate")
  legcol <- c(colpt)
  legpch <- c(pchpt)
  leglwd <- c(lwdpt)
  if (is.element("fkb", inputlist)) {
    plotbp <- TRUE
    legtext <- c("Point estimate", "Bootstrapped estimates")
    legcol <- c(colpt, coli)
    legpch <- c(pchpt, 15)
    leglwd <- c(lwdpt, NA)
    YLIM[2] <- max(x$fkb[ind,])
  }
  if (is.element("fk_ci", inputlist)) {
    plotci <- TRUE
    legtext <- c("Point estimate", paste(100*x$prob, "% confidence interval", sep = ""))
    legcol <- c(colpt, coli)
    legpch <- c(pchpt, NA)
    YLIM[2] <- max(x$fk_ci[,ind])
  }
  plot(k, x$fk[ind], ylim = YLIM, las = las, type = "n", yaxt = "n",
       xlab = "k", ylab = "f(k)", ...)
  tmp <- axTicks(2)
  axis(2, at = tmp[tmp >= 0], las = las)
  abline(h = 0, col = col0, lwd = lwd0)
  if (plotbp) {
    if (sum(ind) == 1) {
      tmp <- x$fkb[ind,]
    } else {
      tmp <- t(x$fkb[ind,])
    }
    boxplot(tmp, at = k, add = TRUE, xaxt = "n", yaxt = "n", boxwex = boxwex,
            col = colibg, border = coli)
    if (plotmu) {
      boxplot(x$mub, at = YLIM[1]/1.5, add = TRUE, xaxt = "n", yaxt = "n", boxwex = boxwex/4,
              col = colibg, border = coli, horizontal = TRUE)
    }
  }
  if (plotci) {
    suppressWarnings(
      arrows(x0 = k, y0 = x$fk_ci[1, ind], y1 = x$fk_ci[2, ind], length = length,
             code = 3, col = coli, angle = 90, lwd = lwdpt)
    )
    if (plotmu) {
      arrows(x0 = x$mu_ci[1], y0 = YLIM[1]/1.5, x1 = x$mu_ci[2], length = length,
             code = 3, col = coli, angle = 90, lwd = lwdpt)
    }
  }
  #Plot point estimates
  lines(k, x$fk[ind], type = "o", pch = pchpt, lwd = lwdpt, col = colpt)
  if (plotmu) {
    axis(2, at = YLIM[1]/1.5, labels = expression(mu), las = las, cex.axis = 1.5)
    points(x$mu, YLIM[1]/1.5, pch = pchpt, lwd = lwdpt, col = colpt)
  }
  if (plotlegend) {
    do.call(legend, append(list(legend = legtext, col = legcol, pch = legpch, lwd = leglwd),
                         legendargs))
  }
}
