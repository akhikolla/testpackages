###Function for plotting a time series of data at each location on the visual field
#'
#' PlotVfTimeSeries
#'
#' Plots a time series at each location of the Humphrey Field Analyzer-II visual field .
#'
#' @param Y a time series variable to be plotted.
#'
#' @param Location a variable corresponding to the location on the visual field
#'  that the time series variable was observed.
#'
#' @param Time a variable corresponding to the time that the time series variable
#'  was observed.
#'
#' @param main an overall title for the plot.
#'
#' @param xlab a title for the x axis.
#'
#' @param ylab a title for the y axis.
#'
#' @param line.col color for the regression line, either character string corresponding
#'  to a color or a integer (default = "red").
#'
#' @param line.reg logical, determines if there are regression lines printed (default = TRUE).
#'
#' @param line.type integer, specifies the type of regression line printed (default = 1).
#'
#' @details \code{PlotVfTimeSeries} is used in the application of glaucoma progression.
#'  In each cell is the observed DLS at each location over visits, with the red line
#'  representing a linear regression trend.
#'
#' @examples
#' data(VFSeries)
#' PlotVfTimeSeries(Y = VFSeries$DLS,
#'                   Location = VFSeries$Location,
#'                   Time = VFSeries$Time,
#'                   main = "Visual field sensitivity time series \n at each location",
#'                   xlab = "Days from baseline visit",
#'                   ylab = "Differential light sensitivity (dB)")
#'
#'
#' @author Samuel I. Berchuck
#'
#' @export
PlotVfTimeSeries <- function(Y, Location, Time,
                             main = "Visual field sensitivity time series \n at each location",
                             xlab = "Time from first visit (days)",
                             ylab = "Sensitivity (dB)",
                             line.col = "red",
                             line.reg = TRUE,
                             line.type = 1) {

  ###Logical function to check for colors
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)),
               error = function(e) FALSE)
    })
  }

  ###Check Inputs
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if (missing(Y)) stop('"Y" is missing')
  if (missing(Location)) stop('"Location" is missing')
  if (missing(Time)) stop('"Time" is missing')
  if (!is.character(main)) stop('"main" must be a character string')
  if (!is.character(xlab)) stop('"xlab" must be a character string')
  if (!is.character(ylab)) stop('"ylab" must be a character string')
  if (!any(areColors(line.col))) stop('"line.col" can only include colors')
  if (!is.logical(line.reg)) stop('"line.reg" must be a logical')
  if (!is.wholenumber(line.type)) stop('"line.type" must be an integer')

  ###Function inputs
  # Y <- YObserved
  # Location <- Location
  # Time <- Time
  # main = "Visual field sensitivity time series \n at each location"
  # xlab = "Time from first visit (days)"
  # ylab = "Sensitivity (dB)"
  pardefault <- suppressWarnings(par(no.readonly = T))

  ###Collect and sort data
  VF <- data.frame(cbind(Location, Time, Y))
  VF <- VF[order(VF$Time), ]
  VF <- VF[order(VF$Location), ]

  ###Compute Summary Statistics
  max.VF <- max(abs(range(Y)))
  max.Time <- max(abs(Time))
  y_breaks <- round(seq(0, 40, by = 10))
  x_breaks <- round(seq(0, 100*(max.Time%/%100 + as.logical(max.Time%%100)), by = 100)) #Round up to the nearest 100th

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
  all <- 1 : max(VF$Location)
  blind_spot <- c(26, 35)
  remaining <- all[-blind_spot]

  ###Plot Time Series at Each Location
  par(mar = c(0, 0, 0, 0), oma = c(5, 10, 10, 5), mgp = c(3, 1, 0))
  for (i in 1 : max(VF$Location)) {
    if (i %in% remaining) {
      ph <- VF[VF$Location == i, ]
      plot(ph[ , 2], ph[ , 3], type = "l", xaxt = "n", yaxt = "n", xlim = c(0, max.Time), ylim = c(0, 40))
      points(ph[ , 2], ph[ , 3], pch = ".")
      if (line.reg) abline(lm(ph[ , 3] ~ ph[ , 2]), col = line.col, lty = line.type)
    }
    if (i %in% blind_spot) plot(ph[ , 2], ph[ , 3], type = "n", xaxt = "n", yaxt = "n", xlim = c(0, max.Time), ylim = c(0, 40))
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
