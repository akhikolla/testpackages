###Function for plotting a time series of data at each location on the visual field
#'
#' PlotCP
#'
#' Plots estimated visual field sensitivities using change point model.
#'
#' @param object a spCP regression object.
#'
#' @param data a dataframe containing the raw sensitivities.
#'
#' @param location a character string indicating the column name in data of the variable of locations on the visual field.
#'
#' @param time a character string indicating the column name in data of the variable of visual field testing times.
#'
#' @param dls a character string indicating the column name in data of the variable of the raw visual field sensitivities.
#'
#' @param main an overall title for the plot.
#'
#' @param xlab a title for the x axis.
#'
#' @param ylab a title for the y axis.
#'
#' @param col color for the regression line, either character string corresponding
#'  to a color or a integer (default = "red").
#'
#' @param ci.col color for the confidence intervals for the regression line, either character string corresponding
#'  to a color or a integer (default = "red").
#'
#' @param cp.col color for the optional change point, either character string corresponding
#'  to a color or a integer (default = "blue").
#'
#' @param cp.ci.col color for the confidence intervals of the optional change point, either character string corresponding
#'  to a color or a integer (default = "blue").
#'
#' @param line logical, determines if there are regression lines printed (default = TRUE).
#'
#' @param ci logical, determines if there are confidence intervals printed (default = TRUE).
#'
#' @param cp.line logical, determines if there the change point is printed (default = FALSE).
#'
#' @param cp.ci logical, determines if there the change point confidence intervals are printed (default = FALSE).
#'
#' @param lwd integer, specifies the width of the regression line (default = 1).
#'
#' @param ci.lwd integer, specifies the width of the confidence intervals (default = 1).
#'
#' @param cp.lwd integer, specifies the width of the change point line (default = 1).
#'
#' @param cp.ci.lwd integer, specifies the width of the change point confidence intervals (default = 1).
#'
#' @param lty integer, specifies the type of regression line (default = 1).
#'
#' @param ci.lty integer, specifies the type of confidence intervals (default = 2).
#'
#' @param cp.lty integer, specifies the type of change point line (default = 1).
#'
#' @param cp.ci.lty integer, specifies the type of change point confidence intervals (default = 2).
#'
#' @details \code{PlotCP} is used in the application of glaucoma progression.
#'  The function is capable of plotting the observed DLS values across the visual field,
#'  along with the estimated mean process (with 95 percent credible intervals) and the
#'  estimated mean posterior change point location (with 95 percent credible intervals).
#'
#' @author Samuel I. Berchuck
#'
#' @export
PlotCP <- function(object,
				   data,
				   location = "Location",
				   time = "Time",
				   dls = "DLS",
 				   line = TRUE,
				   ci = TRUE,
				   lwd = 1,
				   lty = 1,
				   col = 2,
				   ci.lwd = 1,
				   ci.lty = 2,
				   ci.col = 2,
				   cp.line = FALSE,
				   cp.ci = FALSE,
				   cp.lwd = 1,
				   cp.lty = 1,
				   cp.col = 4,
				   cp.ci.lwd = 1,
				   cp.ci.lty = 2,
				   cp.ci.col = 4,
				   main = "Estimated visual field sensitivity\nusing change points",
				   xlab = "Years from baseline visit",
				   ylab = "Sensitivity (dB)") {

  ###Logical function to check for colors
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)),
               error = function(e) FALSE)
    })
  }

  ###Check Inputs
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if (missing(object)) stop('"object" is missing')
  if (!(is.spCP(object))) stop('"object" must be of class spCP')
  if (missing(data)) stop('"data" is missing')
  if (!(is.data.frame(data))) stop('"data" must be of class data.frame')
  if (!is.character(dls)) stop('The variable "dls" must be a character string')
  if (!is.character(location)) stop('The variable "location" must be a character string')
  if (!is.character(time)) stop('The variable "time" must be a character string')
  if (!(dls %in% colnames(data))) stop('The variable name "dls" must be contained in data')
  if (!(location %in% colnames(data))) stop('The variable name "location" must be contained in data')
  if (!(time %in% colnames(data))) stop('The variable name "time" must be contained in data')
  if (!is.character(main)) stop('"main" must be a character string')
  if (!is.character(xlab)) stop('"xlab" must be a character string')
  if (!is.character(ylab)) stop('"ylab" must be a character string')
  if (!any(areColors(col))) stop('"col" can only include colors')
  if (!any(areColors(ci.col))) stop('"ci.col" can only include colors')
  if (!any(areColors(cp.col))) stop('"cp.col" can only include colors')
  if (!any(areColors(cp.ci.col))) stop('"cp.ci.col" can only include colors')
  if (!is.logical(line)) stop('"line" must be a logical')
  if (!is.logical(ci)) stop('"ci" must be a logical')
  if (!is.logical(cp.line)) stop('"cp.line" must be a logical')
  if (!is.logical(cp.ci)) stop('"cp.ci" must be a logical')
  if (!is.wholenumber(lwd)) stop('"lwd" must be an integer')
  if (!is.wholenumber(ci.lwd)) stop('"ci.lwd" must be an integer')
  if (!is.wholenumber(cp.lwd)) stop('"cp.lwd" must be an integer')
  if (!is.wholenumber(cp.ci.lwd)) stop('"cp.ci.lwd" must be an integer')
  if (!is.wholenumber(lty)) stop('"lty" must be an integer')
  if (!is.wholenumber(ci.lty)) stop('"ci.lty" must be an integer')
  if (!is.wholenumber(cp.lty)) stop('"cp.lty" must be an integer')
  if (!is.wholenumber(cp.ci.lty)) stop('"cp.ci.lty" must be an integer')

  ###Save original options
  pardefault <- suppressWarnings(par(no.readonly = T))

  ###Set data variables
  Location <- data[, location]
  Time <- data[, time]
  DLS <- data[, dls]

  ###Read in raw output
  beta0 <- (object$beta0)
  beta1 <- (object$beta1)
  theta <- (object$theta)
  ScaleY <- object$datobj$ScaleY

  ###Number of samples
  NKeep <- dim(theta)[1]
  Day <- unique(Time)
  Nu <- length(Day)

  ###Compute Summary Statistics
  max.VF <- max(abs(range(DLS)))
  max.Time <- max(abs(Time))
  y_breaks <- round(seq(0, 40, by = 10))
  if (max.Time < 7) x_breaks <- round(seq(0, 1*(max.Time%/%1 + as.logical(max.Time%%1)), by = 1)) #Round up to the nearest 1th
  if ((max.Time > 6) & (max.Time < 26)) x_breaks <- round(seq(0, 5*(max.Time%/%5 + as.logical(max.Time%%5)), by = 5)) #Round up to the nearest 5th
  if ((max.Time > 25) & (max.Time < 71)) x_breaks <- round(seq(0, 10*(max.Time%/%10 + as.logical(max.Time%%10)), by = 10)) #Round up to the nearest 10th
  if ((max.Time > 70) & (max.Time < 351)) x_breaks <- round(seq(0, 50*(max.Time%/%50 + as.logical(max.Time%%50)), by = 50)) #Round up to the nearest 50th
  if (max.Time > 350) x_breaks <- round(seq(0, 100*(max.Time%/%100 + as.logical(max.Time%%100)), by = 100)) #Round up to the nearest 100th

  ###Create layout matrix
  layout.matrix<-matrix(c(0,0,0,1,2,3,4,0,0,
                          0,0,5,6,7,8,9,10,0,
                          0,11,12,13,14,15,16,17,18,
                          19,20,21,22,23,24,25,26,27,
                          28,29,30,31,32,33,34,35,36,
                          0,37,38,39,40,41,42,43,44,
                          0,0,45,46,47,48,49,50,0,
                          0,0,0,51,52,53,54,0,0),nrow=8,ncol=9,byrow=TRUE)
  pp<-layout(layout.matrix,rep(1,3),rep(1,9),TRUE)
  # layout.show(pp)

  ###Clarify Blind Spot
  all <- 1 : max(Location)
  blind_spot <- c(26, 35)
  remaining <- all[-blind_spot]
  indeces <- c(1:25, NA, 27:34 - 1, NA, 36:54 - 2)

  ###Plot Time Series at Each Location
  par(mar = c(0, 0, 0, 0), oma = c(5, 10, 10, 5), mgp = c(3, 1, 0))
  for (i in 1 : max(Location)) {
    if (i %in% remaining) {
      plot(-100, -100, type = "l", xaxt = "n", yaxt = "n", xlim = c(0, max.Time), ylim = c(0, 40))
      points(Time[Location == i], DLS[Location == i], pch = 16)
      if (line | ci) {
        Estimate <- matrix(nrow = NKeep, ncol = Nu)
        for (s in 1:NKeep) {
          CP <- theta[s, indeces[i]]
          Estimate[s, ] <- beta0[s, indeces[i]] + beta1[s, indeces[i]] * (Day - CP) * (Day > CP)
        }
        if (line) lines(Day, pmax(0, apply(Estimate, 2, mean)) * ScaleY, col = col, lwd = lwd, lty = lty)
        if (ci) {
          lines(Day, pmax(0, apply(Estimate, 2, function(x) quantile(x, probs = c(0.025)))) * ScaleY, lty = ci.lty, col = ci.col, lwd = ci.lwd)
          lines(Day, pmax(0, apply(Estimate, 2, function(x) quantile(x, probs = c(0.975)))) * ScaleY, lty = ci.lty, col = ci.col, lwd = ci.lwd)
        }
      }
      if (cp.line) abline(v = mean(theta[, indeces[i]]), col = cp.col, lty = cp.lty, lwd = cp.lwd)
      if (cp.ci) {
      	abline(v = quantile(theta[, indeces[i]], prob = 0.025), col = cp.ci.col, lty = cp.ci.lty, lwd = cp.ci.lwd)
		abline(v = quantile(theta[, indeces[i]], probs = 0.975), col = cp.ci.col, lty = cp.ci.lty, lwd = cp.ci.lwd)
      }
    }
    if (i %in% blind_spot) plot(-100, -100, type = "n", xaxt = "n", yaxt = "n", xlim = c(0, max.Time), ylim = c(0, 40))
    if (i %in% c(52, 54)) axis(1, at = x_breaks)
    if (i %in% c(1, 3)) axis(3, at = x_breaks)
    if (i %in% c(5, 19, 37, 51)) axis(2, at = y_breaks, las = 2)
    if (i %in% c(4, 18, 36, 50)) axis(4, at = y_breaks, las = 2)
  }

  ###Add Title
  title(main = list(main, cex = 2.5, col = "black", font = 2),
        xlab = list(xlab, cex = 2, col = "black", font = 1),
        ylab = list(ylab, cex = 2, col = "black", font = 1), outer = TRUE)

  ###Return par to default
  suppressMessages(par(pardefault))

###End function
}
