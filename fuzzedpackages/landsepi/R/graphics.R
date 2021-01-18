# Part of the landsepi R package.
# Copyright (C) 2017 Loup Rimbaud <loup.rimbaud@inrae.fr>
#                    Julien Papaix <julien.papaix@inrae.fr>
#                    Jean-François Rey <jean-francois.rey@inrae.fr>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,i
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#



#' @title Plotting the landscape
#' @name plotland
#' @description Plots a landscape with colors or hatched lines to represent different types of fields
#' @param landscape a spatialpolygon object containing field coordinates
#' @param COL vector containing the color of each field
#' @param DENS vector containing the density of hatched lines for each field
#' @param ANGLE vector containing the angle of hatched lines for each field
#' @param COL.LEG vector containing the colors in the first legend
#' @param DENS.LEG vector containing the density of hatched lines in the second legend
#' @param ANGLE.LEG vector containing the angle of hatched lines in the second legend
#' @param TITLE title of the graphic
#' @param SUBTITLE subtitle of the graphic
#' @param LEGEND1 labels in the first legend (colors)
#' @param LEGEND2 labels in the second legend (hatched lines)
#' @param TITLE.LEG2 title for the second legend
#' @examples
#' \dontrun{
#' ## Draw a landscape with various colours
#' landscapeTEST1
#' plotland(landscapeTEST1,
#'   COL = 1:length(landscapeTEST1),
#'   DENS = rep(0, length(landscapeTEST1)), ANGLE = rep(30, length(landscapeTEST1))
#' )
#' }
#' @importFrom splancs polymap
#' @include RcppExports.R Math-Functions.R
# @S3method plot land
#' @export
plotland <- function(landscape, COL = rep(0, length(landscape)), DENS = rep(0, length(landscape)),
                     ANGLE = rep(30, length(landscape)), COL.LEG = unique(COL), DENS.LEG = unique(DENS),
                     ANGLE.LEG = unique(ANGLE), TITLE = "", SUBTITLE = "", LEGEND1 = rep("", length(COL.LEG)),
                     LEGEND2 = rep("", length(COL.LEG)), TITLE.LEG2 = "") {
  
  ## Bounds of landscape
  bounds <- st_bbox(landscape)

  par(cex = 2, xpd = NA, bg = "white", mar = c(5, 4, 4, 2))
  nPoly <- length(landscape)
  plot(0, 0,
    xlim = c(bounds$xmin, bounds$xmax), ylim = c(bounds$ymin, bounds$ymax), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n", type = "n",
    main = TITLE
  ) # Empty graph
  mtext(SUBTITLE[1], side = 3, line = 0, padj = -.5, las = 1, cex = 1.4)
  if (length(SUBTITLE) > 1) {
    mtext(SUBTITLE[2], side = 3, line = 0, padj = 1.5, las = 1, cex = 1.4)
  }
  for (i in 1:nPoly) {
    polymap(landscape@polygons[[i]]@Polygons[[1]]@coords, add = TRUE, col = COL[i], border = "black", lwd = 2.5)
    polymap(landscape@polygons[[i]]@Polygons[[1]]@coords,
      add = TRUE, col = "black", density = DENS[i],
      angle = ANGLE[i], border = NA
    )
  }
  if (LEGEND1[1] != "") {
    if (TITLE.LEG2 == "") {
      legend(bounds$xmax / 2.66, -bounds$ymax / 40, legend = LEGEND1, fill = COL.LEG, bty = "n")
    } else {
      legend(bounds$xmax / 2.66, -bounds$ymax / 40,
        legend = LEGEND1, col = "black", density = 2 * DENS.LEG, angle = ANGLE.LEG,
        bty = "n"
      )
      legend(-bounds$xmax / 5, bounds$ymax, legend = LEGEND2, fill = COL.LEG, bty = "n", title = TITLE.LEG2)
    }
  }
}



#' @title Plotting allocation of croptypes in a landscape
#' @name plot_allocation
#' @description Plots croptype allocation in the landscape at a given year of the simulation
#' @param landscape a SpatialPolygonsDataFrame
#' @param year year to be plotted
#' @param croptype_names croptype names (for legend)
#' @param title title of the graphic
#' @param subtitle subtitle of the graphic
#' @param filename name of the .png file to be generated
#' @return a png file.
#' @seealso \link{plotland}
#' @examples
#' \dontrun{
#' landscape <- landscapeTEST1
#' croptypes <- data.frame(sample.int(3, length(landscape), replace = TRUE))
#' allocation <- SpatialPolygonsDataFrame(landscape, croptypes, match.ID = TRUE)
#' plot_allocation(allocation, 1,
#'   title = "Simulated landscape", subtitle = "Year 1",
#'   filename = paste(getwd(), "/landscape.png", sep = "")
#' )
#' }
#' @importFrom grDevices dev.off graphics.off png tiff
#' @export
plot_allocation <- function(landscape, year, croptype_names = c(), title = "", subtitle = "", filename = "landscape.png") {
  rotation <- as.matrix(landscape@data)

  rotation2 <- matrix(as.numeric(factor(rotation)), ncol = ncol(rotation)) ## Transform rotation to have croptypes indices as 1, 2, 3, ...
  Ncroptypes <- max(rotation2)
  colfunc <- colorRampPalette(c("black", "white"))
  croptype_col <- colfunc(Ncroptypes + 1) ## +1 to avoid picking white colour
  croptype_col[1] <- "white" ## first croptype
  density <- rep(0, Ncroptypes)
  angle <- rep(0, Ncroptypes)
  if (length(croptype_names) == 0) {
    croptype_names <- as.character(paste("Croptype", levels(factor(rotation))))
  }

  png(filename = filename, width = 1000, height = 1000)
  plotland(
    landscape, croptype_col[rotation2[, year]], density[rotation2[, year]],
    angle[rotation2[, year]], croptype_col[1:Ncroptypes],
    density, angle, title, subtitle, croptype_names
  )
  dev.off()
}



#' @title Plotting pathotype frequencies
#' @name plot_freqPatho
#' @description Plots in a .tiff file the dynamics of pathotype frequencies with respect to
#' pathogen adaptation to a specific resistance gene.
#' @param name_gene a string specifying the name of the gene under investigation
#' @param Nlevels_aggressiveness number of pathotypes with respect to the gene under investigation
#' @param I_aggrProp a matrix giving the frequency of every pathotype (rows) for every time-step (columns)
#' @param nTS number of simulated time-steps
#' @param Nyears number of simulated cropping seasons
#' @param nTSpY number of time-steps per cropping season
#' @examples
#' \dontrun{
#' freqMatrix <- matrix(0, nrow = 2, ncol = 100)
#' freqMatrix[2, 26:100] <- (26:100) / 100
#' freqMatrix[1, ] <- 1 - freqMatrix[2, ]
#' plot_freqPatho(
#'   index_gene = 1,
#'   Nlevels_aggressiveness = 2,
#'   freqMatrix,
#'   nTS = 100,
#'   Nyears = 10,
#'   nTSpY = 10
#' )
#' }
#' @export
plot_freqPatho <- function(name_gene, Nlevels_aggressiveness, I_aggrProp, nTS, Nyears, nTSpY) {
  COL.grey <- gray(0:150 / 150)
  COL.grey <- COL.grey[length(COL.grey):1]
  TITLE <- paste("Pathotype frequencies relative to", name_gene)
  LABELS <- c("\nNon-adapted", "\nAdapted")
  if (Nlevels_aggressiveness > 2) {
    LABELS <- c("\nMaladapted", rep(NA, Nlevels_aggressiveness - 2), "\nAdapted")
  }

  image(
    x = 1:Nlevels_aggressiveness, y = 1:nTS, z = I_aggrProp, col = COL.grey,
    ylab = "", xlab = "Pathotype", main = TITLE, axes = F,
    zlim = c(0, 1), ylim = c(1, nTS + 1), cex.main = 0.8
  )
  box()
  axis(side = 1, at = 1:Nlevels_aggressiveness, labels = LABELS)
  if (Nyears == 1) {
    axis(2, at = round(seq(1, nTS, length.out = 8)), las = 1)
    title(ylab = "Evolutionnary time (days)")
  } else {
    axis(2,
      at = seq(1, nTS + 1, nTSpY * ((Nyears - 1) %/% 10 + 1)),
      labels = seq(0, Nyears, ((Nyears - 1) %/% 10 + 1)), las = 1
    )
    title(ylab = "Evolutionnary time (years)")
  }
}
