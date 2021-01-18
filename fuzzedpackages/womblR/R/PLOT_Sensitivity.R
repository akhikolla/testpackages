###Function used to plot sensitivity values on the visual field
#'
#' PlotSensitivity
#'
#' Plots a heat map of the differential light sensitivity on the Humphrey Field
#' Analyzer-II visual field.
#'
#' @param Y variable to be plotted on the visual field (e.g. differential light sensitivity).
#' @param main an overall title for the plot.
#' @param zlim the limits used for the legend (default are the minimum and maximum of Y).
#' @param color a vector of character strings representing the color palette.
#' @param bins the number of bins used to refine the color palette for the figure and legend.
#' @param legend logical, indicating whether the legend should be present (default = TRUE).
#' @param legend.lab a label for the legend (default = "DLS (dB)").
#' @param legend.round integer, indicating the digits that the legend labels are rounded to
#'  (default = 0).
#' @param border logical, indicating whether there should be a border around the visual field (default = TRUE).
#' @details \code{PlotSensitivity} is used in the application of glaucoma progression to
#'  plot a variable across the visual field in the form of a heat map.
#' @examples
#' data(VFSeries)
#' PlotSensitivity(Y = VFSeries$DLS[VFSeries$Visit == 1],
#'                   main = "Sensitivity estimate (dB) at each \n location on visual field",
#'                   legend.lab = "DLS (dB)",
#'                   zlim = c(10, 35),
#'                   bins = 250)
#' @author Samuel I. Berchuck
#' @export
PlotSensitivity <- function(Y = Y, main = "Sensitivity Estimate (dB) at each \nlocation on visual field",
                            legend.lab = "DLS (dB)", zlim = c(10, 35), bins = 100, border = TRUE, legend = TRUE,
                            color = c("yellow", "orange", "red"), legend.round = 0) {

  ##Note: Depends on library classInt
  # You need the suggested package for this function
  if (!requireNamespace("classInt", quietly = TRUE)) {
    stop("classInt needed for this function to work. Please install it.",
         call. = FALSE)
  }

  ###Check zlim missing
  if (missing(zlim)) zlim <- c(min(Y), max(Y))

  ###Create Legend Cutoffs
  labs <- levels(cut(zlim, bins))
  labs <- cbind(lower = as.numeric(sub("\\((.+),.*","\\1", labs)), upper = as.numeric(sub("[^,]*,([^]]*)\\]","\\1", labs)))
  legvals <- as.numeric(c(labs[1, 1], labs[ , 2]))
  legvals[1] <- -Inf
  legvals[length(legvals)] <- Inf

  ###Get color specification
  colbr <- colorRampPalette(color)
  colpal <- colbr(bins)

  ###Get colors for each observation
  # cuts <- as.character(apply(matrix(Y[!is.na(Y)], ncol = 1), 1, cut, legvals, labels = colpal))
  cuts <- cut(Y[!is.na(Y)], breaks = legvals)
  cuts <- colpal[as.numeric(cuts)]

  ###Create plotting functions
  square <- function(x, y, col) symbols(x, y, squares = 1, fg = col, bg = col, inches = FALSE, add = TRUE)
  format0 <- function(x, legend.round) format(round(x,legend.round),nsmall=legend.round)

  ###Get square coordinates
  Loc <- data.frame(x = c(4:7, 3:8, 2:9, 1:9, 1:9, 2:9, 3:8, 4:7), y = c(rep(1, 4), rep(2, 6), rep(3, 8), rep(4, 9), rep(5, 9), rep(6, 8), rep(7, 6), rep(8, 4)))
  Loc <- Loc[order(Loc$y, decreasing = TRUE),]
  rownames(Loc) <- 1 : 54
  Loc <- Loc[-c(26, 35), ] #remove blind spot

  ###Initiate figure with squares
  pardefault <- suppressWarnings(par(no.readonly = T))
  par(mfcol = c(1, 1), pty = "m", mai = c(0, 0, 0.75, 0))
  # plot(1, 1, main = main, type = "n", yaxt = "n", xaxt = "n", bty = "n", xlim = c(-2, 14), ylim = c(2, 7), asp = 1, ylab = "", xlab = "")
  plot(1, 1, type = "n", yaxt = "n", xaxt = "n", bty = "n", xlim = c(0.5, 13), ylim = c(2, 7), asp = 1, ylab = "", xlab = "")
  title(main = main, cex.main = 1.7)
  for (i in 1 : 52) {
    x <- Loc[i, 1] + 0.5
    y <- Loc[i ,2] + 0.5
    square(x, y, col = cuts[i])
  }
  square(8 + 0.5, 5 + 0.5, col = "grey")
  square(8 + 0.5, 4 + 0.5, col = "grey")

  ###Add border
  if (border) {
    hloop<-list(4:7,c(3,8),c(2,9),1,NULL,1,c(2,9),c(3,8),4:7)
    vloop<-list(4:5,c(3,6),c(2,7),c(1,8),NULL,NULL,NULL,c(1,8),c(2,7),3:6)
    for (j in 1:9) {
      for (i in hloop[[j]]) {
        segments(i,j,i+1,j,lwd = 1.5)
      }
    }
    for (i in 1:10) {
      for (j in vloop[[i]]) {
        segments(i,j,i,j+1,lwd = 1.5)
      }
    }
  }

  ###Add legend
  if (legend) {
    if (missing(zlim)) zlim <- c(min(Y), max(Y))
    NColors <- length(colpal)
    Vertical <- seq(3, 7, length.out = NColors)
    for (i in 1 : NColors) segments(11, Vertical[i], 11.75, Vertical[i], col = colpal[i], lwd = 1.5)
    minx <- zlim[1]
    maxx <- zlim[2]
    LegendPV <- seq(minx, maxx, length.out = 5)
    segments(11.75, 3, 11.75, 7, lwd = 1.5)
    segments(11 ,3 ,11 ,7 , lwd = 1.5)
    segments(11 ,7 ,11.75, 7, lwd = 1.5)
    segments(11 ,3 ,11.75, 3, lwd = 1.5)
    for (i in 1 : length(LegendPV)) {
      text(12.75, (3:7)[i], format0(LegendPV[i], legend.round))
      segments(11.75, (3:7)[i], 12, (3:7)[i], lwd = 1.5)
    }
    text(11.5, 7.5, legend.lab)
  }

  ###Return to default par setting
  par(pardefault)

###End function
}
