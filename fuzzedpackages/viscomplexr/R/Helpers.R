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



#' Adjust ylim to xlim
#'
#' This simple function is useful for adjusting x and y coordinate ranges
#' \code{xlim} and \code{ylim} in order to maintain a desired display ratio. The
#' former must be given, the latter will be adjusted.
#'
#' For certain purposes, e.g. producing a graph that exactly matches a screen,
#' the x and y coordinates must be adjusted to match a given display ratio. If
#' the horizontal range, \code{xlim}, the desired ratio, \code{x_to_y} and the
#' desired center of the y-range, \code{centerY} are provided, this function
#' returns an adapted vertical range, that can be used as \code{ylim} in any
#' plot including \code{\link{phasePortrait}}.
#'
#' @param xlim Numeric vector of length 2; the fixed lower and upper boundary
#'   of the horizontal coordinate range
#'
#' @param centerY The vertical coordinate which the output range is to be
#'   centered around (default = 0)
#'
#' @param x_to_y The desired ratio of the horizontal (x) to the vertical (y)
#'   range. Default is 16/9, a display ratio frequently used for computer or
#'   mobile screens
#'
#' @return A numeric vector of length 2; the lower and upper boundary of the
#'   resulting vertical coordinate range
#'
#' @family helpers
#'
#' @export
#'
#' @examples
#' # Make a phase portrait of a Jacobi theta function that fully covers a
#' # plot with a display aspect ratio of 4/3.
#' # 10 inch wide window with 4/3 display ratio (x/y)
#' \donttest{
#' # x11(width = 10, height = 10 * 3/4) # Screen device commented out
#'                                      # due to CRAN test requirements.
#'                                      # Use it when trying this example
#' xlim <- c(-3, 3)
#' ylim <- ylimFromXlim(xlim, centerY = -0.3, x_to_y = 4/3)
#' op <- par(mar = c(0, 0, 0, 0), bg = "black") # Omit all plot margins
#' phasePortrait(jacobiTheta, moreArgs = list(tau = 1i/2 - 1/3),
#'   xlim = xlim, ylim = ylim, # Apply the coordinate ranges
#'   xaxs = "i", yaxs = "i",   # Allow for now room between plot and axes
#'   nCores = 1) # Max. two cores allowed on CRAN
#'               # not a limit for your own use
#' par(op)
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
ylimFromXlim <- function(xlim, centerY = 0, x_to_y = 16/9) {

  yRange <- abs(diff(xlim)) / x_to_y
  ylim   <- centerY + c(-1/2, 1/2) * yRange

  return(ylim)
}



#' Adjust xlim to ylim
#'
#' This simple function is useful for adjusting x and y coordinate ranges
#' \code{xlim} and \code{ylim} in order to maintain a desired display ratio. The
#' latter must be given, the former will be adjusted.
#'
#' For certain purposes, e.g. producing a graph that exactly matches a screen,
#' the x and y coordinates must be adjusted to match a given display ratio. If
#' the vertical range, \code{ylim}, the desired ratio, \code{x_to_y} and the
#' desired center of the x-range, \code{centerX}, are provided, this function
#' returns an adpated vertical range, that can be used as \code{ylim} in any
#' plot including \code{\link{phasePortrait}}.
#'
#' @param ylim Numeric vector of length 2; the fixed lower and upper boundary
#'   of the vertical coordinate range
#'
#' @param centerX The horizontal coordinate which the output range is to be
#'   centered around (default = 0)
#'
#' @param x_to_y The desired ratio of the horizontal (x) to the vertical (y)
#'   range. Default is 16/9, a display ratio frequently used for computer or
#'   mobile screens
#'
#' @return A numeric vector of length 2; the lower and upper boundary of the
#'   resulting vertical coordinate range
#'
#' @family helpers
#'
#' @export
#'
#' @examples
#' # Make a phase portrait of a pretty function that fully covers a
#' # plot with a display aspect ratio of 5/4.
#'
#' # 9 inch wide window with 5/4 display ratio (x/y)
#' \donttest{
#' # x11(width = 9, height = 9 * 4/5) # Screen device commented out
#'                                    # due to CRAN test requirements.
#'                                    # Use it when trying this example
#' ylim <- c(-8, 7)
#' xlim <- xlimFromYlim(ylim, centerX = 0, x_to_y = 5/4)
#' op <- par(mar = c(0, 0, 0, 0), bg = "black") # Omit all plot margins
#' phasePortrait("exp(cosh(1/(z - 2i + 2)^2 * (1/2i - 1/4 + z)^3))", pType = "pm",
#' xlim = xlim, ylim = ylim, # Apply the coordinate ranges
#' xaxs = "i", yaxs = "i",   # Allow for now room between plot and axes
#' nCores = 2) # Max. two cores allowed on CRAN
#'             # not a limit for your own use
#' par(op)
#'   \dontshow{
#'   # R CMD check: make sure any open connections are closed afterward
#'   foreach::registerDoSEQ()
#'   doParallel::stopImplicitCluster()
#'   }
#' }
#'
#'
xlimFromYlim <- function(ylim, centerX = 0, x_to_y = 16/9) {

  xRange <- abs(diff(ylim)) * x_to_y
  xlim   <- centerX + c(-1/2, 1/2) * xRange

  return(xlim)
}



