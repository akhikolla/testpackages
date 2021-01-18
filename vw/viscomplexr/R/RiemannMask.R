# R package viscomplexr - phase portraits of functions in the
# complex number plane
# Copyright (C) 2020 Peter Biber
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>



#' Plot a Riemann sphere mask over a phase portrait
#'
#' The function \code{riemannMask} can be used for laying a circular mask over
#' an existing phasePortrait (as generated with the function
#' \code{\link{phasePortrait}}). This mask shades the plot region outside the
#' unit circle. The unshaded area is a projection on the southern or northern
#' Riemann hemisphere. The standard projection used by
#' \code{\link{phasePortrait}}, i.e. \code{invertFlip = FALSE} hereby
#' corresponds to the southern Riemann hemisphere with the origin being the
#' south pole. If \code{\link{phasePortrait}} was called with \code{invertFlip =
#' TRUE}, then the unit circle contains the northern Riemann hemisphere with the
#' point at infinity in the center (see the vignette for more details). Options
#' for adding annotation, landmark points  are available
#' (\insertCite{@see @wegert_visualcpx_2012;textual}{viscomplexr}, p. 41).
#' Several parameters are on hand for adjusting the mask's transparency, color,
#' and similar features. some details, this function behaves less nicely under
#' Windows than under Linux (see Details).
#'
#' There is, unfortunately, a somewhat different behavior of this function under
#' Linux and Windows systems. Under Windows, the region outside the unit circle
#' is only shaded if the whole unit circle fits into the plot region. If only a
#' part of the unit circle is to be displayed, the shading is completely omitted
#' under Windows (annotation etc. works correctly, however), while it works
#' properly on Linux systems. Obviously, the function \code{\link{polypath}},
#' which we are using for creating the unit circle template, is interpreted
#' differently on both systems.
#'
#' @param colMask Color for the shaded area outside the unit circle. Defaults to
#'   "white". Can be any kind of color definition R accepts. I recommend,
#'   however, to use a color definition without a transparency value, because
#'   this would be overridden by the parameter \code{alphaMask}.
#'
#' @param alphaMask Transparency value for the color defined with
#'   \code{colMask}. Has to be a value between 0 (fully transparent) and 1
#'   (totally opaque). Defaults to 0.5.
#'
#' @param circOutline Boolean - if \code{TRUE}, the outline of the unit circle
#'   is drawn. Defaults to
#'   \code{TRUE}.
#'
#' @param circLwd Line width of the unit circle outline. Obviously relevant
#'   only when \code{circOutline == TRUE}. Defaults to 1.
#'
#' @param circleSteps Number of vertices to draw the circle. Defaults to 360
#'   (one degree between two vertices).
#'
#' @param circleCol Color of the unit circle, default is the default foreground
#'   color (\code{par("fg")}).
#'
#' @param gridCross Boolean - if \code{TRUE}, a horizontal and a vertical gray
#'   line will be drawn over the plot region, intersection in the center of the
#'   unit circle. Defaults to \code{FALSE}.
#'
#' @param annotSouth Boolean - add landmark points and annotation for a
#'   \emph{southern} Riemann hemisphere, defaults to \code{FALSE}. This
#'   annotation fits to an image that has been created with
#'   \code{\link{phasePortrait}} and the option \code{invertFlip = FALSE}.
#'
#' @param annotNorth Boolean - add landmark points and annotation for a
#'   \emph{northern} Riemann hemisphere, defaults to \code{FALSE}. This
#'   annotation fits to an image that has been created with
#'   \code{\link{phasePortrait}} and the option \code{invertFlip = TRUE}.
#'
#' @param xlim,ylim optional, if provided must by numeric vectors of length 2
#'   defining plot limits as usual. They define the outer rectangle of the
#'   Riemann mask. If \code{xlim} or \code{ylim} is not provided (the standard
#'   case), the coordinates of the plot window as given by \code{par("usr")}
#'   will be used for the missing component.
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#' # Tangent with fully annotated Riemann masks.
#' # The axis tick marks on the second diagram (Northern hemisphere)
#' # have to be interpreted as the real and imaginary parts of 1/z
#' # (see vignette). The axis labels in this example have been adapted
#' # accordingly.
#' \donttest{
#' # x11(width = 16, height = 8) # Screen device commented out
#'                               # due to CRAN test requirements.
#'                               # Use it when trying this example
#' op <- par(mfrow = c(1, 2), mar = c(4.7, 4.7, 3.5, 3.5))
#' phasePortrait("tan(z)", pType = "pma",
#'               main = "Southern Riemann Hemisphere",
#'               xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
#'               xlab = "real", ylab = "imaginary",
#'               xaxs = "i", yaxs = "i",
#'               nCores = 2) # Max. two cores on CRAN, not a limit for your use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' riemannMask(annotSouth = TRUE, gridCross = TRUE)
#'
#' phasePortrait("tan(z)", pType = "pma",
#'               main = "Northern Riemann Hemisphere",
#'               invertFlip = TRUE,
#'               xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
#'               xlab = "real (1/z)", ylab = "imaginary (1/z)",
#'               xaxs = "i", yaxs = "i",
#'               nCores = 2) # Max. two cores on CRAN, not a limit for your use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' riemannMask(annotNorth = TRUE, gridCross = TRUE)
#' par(op)
#' }
#'
#' # Rational function with Riemann masks without annotation.
#' # The axis tick marks on the second diagram (Northern hemisphere)
#' # have to be interpreted as the real and imaginary parts of 1/z
#' # (see vignette). The axis labels in this example have been adapted
#' # accordingly.
#' \donttest{
#' # x11(width = 16, height = 8) # Screen device commented out
#'                               # due to CRAN test requirements.
#'                               # Use it when trying this example
#' op <- par(mfrow = c(1, 2), mar = c(4.7, 4.7, 3.5, 3.5))
#' phasePortrait("(-z^17 - z^15 - z^9 - z^7 - z^2 - z + 1)/(1i*z - 1)",
#'               pType = "pma",
#'               main = "Southern Riemann Hemisphere",
#'               xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
#'               xlab = "real", ylab = "imaginary",
#'               xaxs = "i", yaxs = "i",
#'               nCores = 2) # Max. two cores on CRAN, not a limit for your use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' riemannMask(annotSouth = FALSE, gridCross = FALSE, circOutline = FALSE,
#'             alphaMask = 0.7)
#'
#' phasePortrait("(-z^17 - z^15 - z^9 - z^7 - z^2 - z + 1)/(1i*z - 1)",
#'               pType = "pma",
#'               main = "Northern Riemann Hemisphere",
#'               invertFlip = TRUE,
#'               xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
#'               xlab = "real (1/z)", ylab = "imaginary (1/z)",
#'               xaxs = "i", yaxs = "i",
#'               nCores = 2) # Max. two cores on CRAN, not a limit for your use
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' riemannMask(annotNorth = FALSE, gridCross = FALSE, circOutline = FALSE,
#'             alphaMask = 0.7)
#' par(op)
#' }
#'
#'
#' @export
#'
#'
riemannMask <- function(colMask     = "white",
                        alphaMask   = 0.5,
                        circOutline = TRUE,
                        circLwd     = 1,
                        circleSteps = 360,
                        circleCol   = par("fg"),
                        gridCross   = FALSE,
                        annotSouth  = FALSE,
                        annotNorth  = FALSE,
                        xlim        = NULL,
                        ylim        = NULL) {

  # Get user plot coordinate extremes and calculate
  # the widths in both directions (x and y)
  coord  <- par("usr")
  xlmt   <- c(coord[1],  coord[2])
  ylmt   <- c(coord[3],  coord[4])

  # Overwrite if xlim and/or ylim are not NULL
  if(!is.null(xlim)) xlmt <- xlim
  if(!is.null(ylim)) ylmt <- ylim

  # Define the outer frame of the mask
  frame  <- list(x = c(xlmt[1], xlmt[1], xlmt[2], xlmt[2]),
                 y = c(ylmt[1], ylmt[2], ylmt[2], ylmt[1]))

  # Define the inner unit circle
  innerCircle <- list(x = cos(seq(360/circleSteps, 360,
                                  by = 360/circleSteps)*2*pi/360),
                      y = sin(seq(360/circleSteps, 360,
                                  by = 360/circleSteps)*2*pi/360))

  # Plot it
  polypath(x = c(frame$x, NA, innerCircle$x),
           y = c(frame$y, NA, innerCircle$y),
           col  = scales::alpha(colMask, alpha = alphaMask),
           rule = "evenodd",
           border = NA)

  # Add circle
  if(circOutline) {
    plotrix::draw.circle(0, 0, 1, nv = circleSteps,
                         lwd = circLwd, border = circleCol)
  }

  # Add grid cross if desired
  if(gridCross) {
    abline(h = 0, col = "grey")
    abline(v = 0, col = "grey")
  }

  # Add annotation if desired (for southern hemisphere)
  if(annotSouth) {
    pts <- list(x = c(-1, 0, 0, 0, 1), y = c(0, 1, 0, -1, 0))
    points(pts$x, pts$y, pch = 21, col = "black",
           bg = c("white", "white", "black", "white", "white"))

    labels <- c("-1", "i", "0", "-i", "1")
    i <- c(1:5)
    lapply(i, function(i, pts, labels) {
      adx    <- switch(labels[i],
        "-1" =   1.5, "i" =  -2.5, "0" = -2.0, "-i" = -0.6, "1" = -1.8)
      ady    <- switch(labels[i],
        "-1" =  -1.2, "i" =  -1.2, "0" = -1.2, "-i" =  2.2, "1" = -1.2)
      text(pts$x[i], pts$y[i], labels[i], vfont = c("serif", "bold"),
           adj = c(adx, ady))
    }, pts = pts, labels = labels)
  } # if annotSouth

  # Add annotation if desired (for northern hemisphere)
  if(annotNorth) {
    pts <- list(x = c(-1, 0, 0, 0, 1), y = c(0, 1, 0, -1, 0))
    points(pts$x, pts$y, pch = 21, col = "black",
           bg = c("white"))

    labels <- c("-1", "i", " ", "-i", "1")
    i <- c(1:5)
    lapply(i, function(i, pts, labels) {
      adx    <- switch(labels[i],
        "-1" =  -0.6, "i" =  -2.5, " " = -2.0, "-i" = -0.6, "1" =  2.6)
      ady    <- switch(labels[i],
        "-1" =  -1.2, "i" =  -1.2, " " = -1.2, "-i" =  2.2, "1" = -1.2)
      text(pts$x[i], pts$y[i], labels[i], vfont = c("serif", "bold"),
           adj = c(adx, ady))
    }, pts = pts, labels = labels)

    text(0, 0, expression(infinity), adj = c(-0.7, -1.0), cex = 1.7, font = 2)
  } # if annotNorth

} # function riemannMask

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

