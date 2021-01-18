#'@name hist.FDX
#'@title Histogram of Raw p-Values
#'
#'@description
#'Computes a histogram of the raw p-values of a \code{FDX} object.
#'
#'@param x          an object of class "\code{FDX}".
#'@param breaks     as in \code{\link{hist}}; here, the Friedman-Diaconis
#'                  algorithm(\code{"FD"}) is used as default.
#'@param main       main title. If \code{NULL} (default), a description string
#'                  is used.
#'@param xlab,ylab  labels for x and y axis.
#'@param plot       a boolean If \code{TRUE} (the default), a histogram is
#'                  plotted, otherwise a list of breaks and counts is returned.
#'@param ...        further arguments to \code{\link{hist}} or
#'                  \code{\link{plot.histogram}}, respectively.
#'
#'@details
#'If \code{x} contains results of a weighted approach, a histogram of the
#'weighted p-values is constructed. Otherwise, it is constituted by the
#'raw ones. 
#'
#'@return
#'An object of class \code{histogram}.
#'
#'@template example
#'@examples
#'
#'DGR <- DGR(raw.pvalues, pCDFlist)
#'hist(DGR)
#'
#'@method hist FDX
#'@importFrom graphics hist
#'@export
hist.FDX <- function(x, breaks = "FD", main = NULL, xlab = NULL, ylab = NULL, plot = TRUE, ...){
  # check if 'main' or 'xlab' arguments are NULL and set defaults accordingly
  if(!grepl("Weighted", x$Method)){
    if(is.null(main))
      main <- "Histogram of raw p-values"
    if(is.null(xlab))
      xlab <- "Raw p-values"
    y <- x$Data$raw.pvalues
  }else{
    if(is.null(main))
      main <- "Histogram of weighted p-values"
    if(is.null(xlab))
      xlab <- "Weighted p-values"
    y <- x$Weighted
  }
  
  # call 'hist' function with raw p.values (default breaks: "FD"); "..." passes
  # all additional 'hist' arguments to this call
  r <- hist(y, breaks = breaks, main = main, xlab = xlab, plot = plot,  ...)
  
  r$xname <- deparse(substitute(x))
  
  if(plot) return(invisible(r)) else return(r)
}


#'@name plot.FDX
#'@title Plot Method for \code{FDX} objects
#'
#'@description
#'Plots raw p-values of a \code{FDX} object and highlights rejected and
#'accepted p-values. If present, the critical values are plotted, too.
#'
#'@param x          an object of class "\code{FDX}".
#'@param col        a numeric or character vector of length 3 indicating the
#'                  colors of the \enumerate{
#'                    \item rejected p-values
#'                    \item accepted p-values
#'                    \item critical values (if present).
#'                  }
#'@param pch        a numeric or character vector of length 3 indicating the
#'                  point characters of the \enumerate{
#'                    \item rejected p-values
#'                    \item accepted p-values
#'                    \item critical values (if present and \code{type.crit}
#'                          is a plot type like \code{'p'}, \code{'b'} etc.).
#'                  }
#'@param lwd        a numeric vector of length 3 indicating the thickness of
#'                  the points and lines.
#'@param type.crit  1-character string giving the type of plot desired for the
#'                  critical values (e.g.: \code{'p'}, \code{'l'} etc; see
#'                  \code{\link{plot}}).
#'@param legend     if NULL, no legend is plotted; otherwise expecting a
#'                  character string like "topleft" etc. or a numeric vector
#'                  of two elements indicating (x, y) coordinates.
#'@param ...        further arguments to \code{\link{plot.default}}.
#'
#'@details
#'If \code{x} contains results of a weighted approach, the Y-axis of the plot
#'is derived from the weighted p-values. Otherwise, it is constituted by the
#'raw ones. 
#'
#'@template example
#'@examples
#'
#'DLR.sd.fast <- DLR(raw.pvalues, pCDFlist)
#'DLR.sd.crit <- DLR(raw.pvalues, pCDFlist, critical.values = TRUE)
#'DLR.su.fast <- DLR(raw.pvalues, pCDFlist, direction = "su")
#'DLR.su.crit <- DLR(raw.pvalues, pCDFlist, direction = "su", critical.values = TRUE)
#'
#'plot(DLR.su.fast)
#'plot(DLR.su.crit, xlim = c(1, 5), ylim = c(0, 0.4))
#'plot(DLR.sd.fast, col = c(2, 4), pch = c(2, 3), lwd = c(2, 2), 
#'     legend = "topleft", xlim = c(1, 5), ylim = c(0, 0.4))
#'plot(DLR.sd.crit, col = c(2, 4, 1), pch = c(1, 1, 4), lwd = c(1, 1, 2), 
#'     type.crit = 'o', legend = c(1, 0.4), lty = 1, xlim = c(1, 5), 
#'     ylim = c(0, 0.4))
#'
#'@method plot FDX
#'@importFrom graphics plot lines points
#'@importFrom methods hasArg
#'@export
plot.FDX <- function(x, col = c(2, 4, 1), pch = c(1, 1, 1), lwd = c(1, 1, 1), type.crit = 'b', legend = NULL, ...){
  # necessary to appease automated R CMD check on CRAN
  main <- ylab <- lty <- NULL
  if(!grepl("Weighted", x$Method)){
    y <- x$Data$raw.pvalues
  }else{
    y <- x$Weighted
  }
  # determine number of tests and rejections
  m <- length(x$Data$raw.pvalues)
  k <- length(x$Indices)
  
  # replicate shorter plot parameter vectors to avoid errors
  col <- rep_len(col, 3)
  pch <- rep_len(pch, 3)
  lwd <- rep_len(lwd, 3)
  
  # get values of ...-arguments
  lst <- list(...)
  
  # start plotting with empty area
  if(!hasArg(main) && !hasArg(ylab))
    plot(y, col = NA, main = x$Method, ylab = paste("Sorted", ifelse(grepl("Weighted", x$Method), "weighted", "raw"),"p-values"), ...)
  else if(!hasArg(main))
    plot(y, col = NA, main = x$Method, ...)
  else if(!hasArg(ylab))
    plot(y, col = NA, ylab = paste("Sorted", ifelse(grepl("Weighted", x$Method), "weighted", "raw"),"p-values"), ...)
  else plot(y, col = NA, ...)
  
  # plot critical values (if present and plotting is requested by the user)
  if(exists('Critical.values', where = x) && type.crit != 'n'){
    lines(x$Critical.values, col = col[3], lwd = lwd[3], pch = pch[3], type = type.crit, ...)
  }
  # plot rejected p-values
  points(sort(y)[1:k], col = col[1], pch = pch[1], lwd = lwd[1], ...)
  # plot accepted p-values
  points((k + 1):m, sort(y[-x$Indices]), col = col[2], pch = pch[2], lwd = lwd[2], ...)
  
  # plot legend
  if(!is.null(legend)){
    n <- 2 + exists('Critical.values', where = x)
    lt <- rep(0, n)
    if(n > 2 && hasArg(lty) && !(type.crit %in% c('p', 'n'))) lt[3] <- lst$lty
    if(length(legend) == 1){
      legend(legend, NULL, c("Rejected", "Accepted", "Critical values")[1:n], col = col[1:n], pch = pch[1:n], lty = lt[1:n], lwd = lwd[1:n])
    }else if(length(legend) == 2){
      legend(legend[1], legend[2], c("Rejected", "Accepted", "Critical values")[1:n], col = col[1:n], pch = pch[1:n], lty = lt[1:n], lwd = lwd[1:n])
    }else warning("Expecting character string or numeric vector of one or two elements for creating a legend")
  }
}


#'@name rejection.path
#'@title Rejection Path Plot (for \code{FDX} objects)
#'
#'@description
#'Displays the number of rejections of the raw p-values in a \code{FDX}
#'object in dependence of the exceedance probability \code{zeta}.
#'
#'@param x                 an object of class "\code{FDX}".
#'@param xlim              the x limits of the plot. If \code{NULL} (default),
#'                         the (0, 1) range is used.
#'@param ylim              the y limits of the plot. If \code{NULL} (default),
#'                         the double of the median of the number of possible
#'                         rejections is used as upper limit.
#'@param main              main title. If \code{NULL} (default), a description
#'                         string is used.
#'@param xlab,ylab         labels for x and y axis.
#'@param verticals         logical; if \code{TRUE}, draw vertical lines at
#'                         steps.
#'@param pch               jump point character.
#'@param ref.show          logical; if \code{TRUE} a vertical reference line
#'                         is plotted, whose height is the number of
#'                         rejections of the original Benjamini-Hochberg (BH)
#'                         procedure.
#'@param ref.col           color of the reference line.
#'@param ref.lty,ref.lwd   line type and thickness for the reference line.
#'@param ...               further arguments to \code{\link{plot.stepfun}}.
#'
#'@return
#'Invisibly returns a \code{stepfun} object that computes the number of
#'rejectionsin dependence on the exceedance probability \code{zeta}.
#'
#'@template example
#'@examples
#'
#'DLR <- DLR(raw.pvalues, pCDFlist)
#'NDLR <- NDLR(raw.pvalues, pCDFlist)
#'
#'rejection.path(DLR, xlim = c(0, 1), ref.show = TRUE, ref.col = "green", ref.lty = 4)
#'rejection.path(NDLR, col = "red", add = TRUE)
#'
#'@importFrom graphics plot abline mtext
#'@importFrom stats stepfun ecdf p.adjust
#'@export
rejection.path <- function(x, xlim = NULL, ylim = NULL, main = NULL, xlab = expression(zeta), ylab = "Number of Rejections", verticals = FALSE, pch = 19, ref.show = FALSE, ref.col = "gray", ref.lty = 2, ref.lwd = 2, ...){
  if(class(x) != "FDX") stop("'x' must be an object of class 'FDX'!")
  # number of hypotheses
  m <- length(x$Data$raw.pvalues)
  # number of BH rejections
  num.rejections.BH <- sum(p.adjust(x$Data$raw.pvalues, "BH") <= x$FDP.threshold)
  
  # create step function
  ecdf.env <- environment(ecdf(x$Adjusted))
  stepfun.x <- stepfun(ecdf.env$x, c(0, m * ecdf.env$y))
  
  lst <- list(...)
  subt <- NULL
  
  # plot
  if(is.null(xlim)) xlim <- c(0, 1)
  if(is.null(ylim)) ylim <- c(0, min(m, 2 * stepfun.x(0.5)))
  if(is.null(main)){
    main <- bquote(bold("Rejection path for"~alpha==.(as.character(x$FDP.threshold))))
    if(!exists('add', where = lst) || (exists('add', where = lst) && !lst$add)) subt <- x$Method
  }
  
  plot(stepfun.x, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab, verticals = verticals, pch = pch, ...)
  if(!is.null(subt)) mtext(subt, line = 0.3)
  
  if(ref.show) abline(h = num.rejections.BH, col = ref.col, lty = ref.lty, lwd = ref.lwd)
  
  return(invisible(stepfun.x))
}
