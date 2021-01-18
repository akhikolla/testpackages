#' Mean residual life plot
#'
#' Calculate mean residual life and plot it to aid the identification of a
#' threshold over which to fit a generalized Pareto distribution
#'
#' Threshold choice for the fitting of the GPD is guided by the shape of the
#' Mean Residual Life plot.  A threshold which is suitably high will have a
#' corresponding mrl plot which is approximately linear in shape above the
#' threshold (up to sampling variation).
#'
#' @aliases mrl print.mrl summary.mrl print.summary.mrl plot.mrl ggplot.mrl
#' @usage mrl(data, umin = min(data), umax = max(data) - 0.1, nint = 100,
#' alpha=.050)
#' \method{print}{mrl}(x, ...)
#' \method{print}{summary.mrl}(x, ...)
#' \method{summary}{mrl}(object, ...)
#' \method{plot}{mrl}(x, xlab="Threshold", ylab="Mean excess", ...)
#' \method{ggplot}{mrl}(data, mapping, xlab = "Threshold",
#'   ylab = "Mean excess", main=NULL,fill="orange", col="blue",
#'   rug=TRUE, addNexcesses=TRUE, textsize=4, ..., environment)
#' @param data A numeric vector.
#' @param umin The minimum value over which to threshold the data.
#' @param umax The maximum value over which to threshold the data.
#' @param alpha Used to determine coverage of confidence interval to plot.
#' Defaults to plotting a 95\% interval.
#' @param nint The number of points at which to compute the plot.
#' @param x,object Arguments to print, summary and plot functions.
#' @param xlab Label for the x-axis. Defaults to \code{xlab="Threshold"}.
#' @param ylab Label for the y-axis. Defaults to \code{ylab="Mean excess"}.
#' @param \dots Optional arguments to \code{plot}.
#' @param col Colour of the line on the MRL plot.
#' @param rug Whether to add raw data as a rug along axis of plot.
#' @param fill Colour of the pointwise confidence region on the MRL plot.
#' @param main Main title.
#' @param addNexcesses Whether to annotate the plot with the numbers of
#' excesses over increasing thresholds. Defaults to \code{addNexcesses=TRUE}.
#' @param textsize Size of text on the plot (ggplot). Defaults to
#' \code{textsize=4}.
#' @param mapping,environment Not used.
#' @return A list with two components. \code{data} is the original data,
#' \code{mrl} is a matrix containing information to produce the mean residual
#' life plot.
#' @note The function was originally written by Stuart Coles and appears in the
#' \code{ismev} package. This version modified by Harry Southworth to allow
#' more control over the appearance of the plot.
#' @author Janet E. Heffernan, Harry Southworth
#' @references S. Coles, An Introduction to Statistical Modeling of Extreme
#' Values, Springer, 2001
#' @keywords models
#' @export mrl
mrl <-
function (data, umin = min(data), umax = max(data) - 0.1,
          nint = 100, alpha=.050) {
    data <- c(data)
    AllData <- data
    x <- lo <- hi <- numeric(nint)
    Threshold <- seq(umin, umax, length = nint)
    z <- qnorm(1 - alpha/2)
    for (i in 1:nint) {
        data <- data[data > Threshold[i]]
        x[i] <- mean(data - Threshold[i])
        sdev <- sqrt(var(data))
        lo[i] <- x[i] - z*sdev/sqrt(length(data))
        hi[i] <- x[i] + z*sdev/sqrt(length(data))
    }

    res <- cbind(threshold=Threshold, MRL=x, lo=lo, hi=hi)
    res <- list(mrl=res, data=AllData)
    oldClass(res) <- 'mrl'
    res
}

#' @export
print.mrl <- function(x, ...){
    print(x$mrl)
    invisible(x)
}

#' @export
summary.mrl <- function(object, ...){
    obj <- list(table=summary(object$mrl))
    oldClass(obj) <- "summary.mrl"
    obj
}

#' @export
print.summary.mrl <- function(x, ...){
    print(x$table)
    invisible(x)
}

#' @export
plot.mrl <- function(x, xlab="Threshold", ylab="Mean excess", ...){

    data <- x$data
    x <- x$mrl

    th <- x[, "threshold"]
    mrl <- x[, "MRL"]
    xl <- x[, "lo"]
    xu <- x[, "hi"]

    plot(th, mrl, type = "l", xlab = xlab, ylab = ylab,
        ylim = c(min(xl[!is.na(xl)]), max(xu[!is.na(xu)])), ...)
    lines(th[!is.na(xl)], xl[!is.na(xl)], lty = 2)
    lines(th[!is.na(xu)], xu[!is.na(xu)], lty = 2)
    rug(data)
    invisible()
}

