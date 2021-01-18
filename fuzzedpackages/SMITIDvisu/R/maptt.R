# This file is part of SMITIDvisu package.
# Copyright (C) 2018-2019 Jean-François Rey <jean-francois.rey@inra.fr>
#
# SMITIDvisu is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMITIDvisu is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SMITIDvisu. If not, see <https://www.gnu.org/licenses/>.
#


#' maptt
#' @description Display a Transmission Tree over a map.
#' @param data Either a data frame that will be converted to a GeoJSON
#' collection, or a string describing a valid GeoJSON collection.
#' The data frame must contain at least columns 'id', 'time', 'X' and 'Y'.
#' It can contain columns 'infectedby', 'probabilities'.
#' Additionnal columns will be added as properties, but will do nothing
#' in this implementation of maptt.
#' See the `df2geojson` function for more informations.
#' @param multipleValuesByTime Vector of strings indicating the df columns names which can contain several values by time. Typically, you would use `c('infectedby','probabilities')` if you have these values.
#' @param circleRadius Numeric value specifying the radius of the nodes in pixels.
#' @param defaultNodeColor String indicating the default color of nodes, if their status doesn't match whith any color. Colors can be specified in hex.
#' @param nodeColorByState List of strings, indicating the color scheme for each node state.
#' @param moveEdgeColor String indicating the color of the edges representing the move of a node.
#' @param color1 String indicating the color corresponding to the minWeight value.
#' @param color2 String indicating the color corresponding to the maxWeight value.
#' @param nbColors Number of colors for the color scheme using a gradient between color1 and color2. These colors will be used to represent the infection edges according to the infection probability. If no probability is used, the edge will use color2. Three intervals are created :
#' color1 will be used for the probabilities between minWeight and weight1.
#' Colors between color1 and color2 will be used for probabilities between weight1 and weight2.
#' color2 will be used for probabilities between weight2 and maxWeight.
#' This setting can be modified directly on the map if 'gradientControl' is activated.
#' @param minWeight Minimal weight.
#' @param maxWeight Maximal weight.
#' @param weight1 Lowest weight for the color scheme. This setting can be modified directly on the map if 'gradientControl' is activated.
#' @param weight2 Greatest weight for the color scheme. This setting can be modified directly on the map if 'gradientControl' is activated.
#' @param autoFocus Boolean indicating if the map should focus at the displayed features at each time. This setting can be toggled directly on the map if 'optionsControl' is activated.
#' @param keepOldFeatures Boolean indicating if old features should be displayed or not. Features are considered "old" if their last 'time' is prior to the current time displayed. This setting can be toggled directly on the map if 'optionsControl' is activated.
#' @param optionsControl Boolean indicating if the options control should be displayed or not
#' @param gradientControl Boolean indicating if the gradient control should be displayed or not
#' @param legend Boolean indicating if the legend should be displayed or not
#' @param width Numeric width for the area in pixels.
#' @param height Numeric hieght for the area in pixels.
#' @param elementId The element ID where the map is displayed
#' @examples
#' library(SMITIDvisu)
#' data(transmissiontree)
#'
#' maptt(tt.events, multipleValuesByTime = c('infectedby', 'probabilities'))
#'
#' # In this example:
#' # - values lower than 20 will be yellow ;
#' # - values between 20 and 25 will use colors between yellow and red ;
#' # - values greater than 25 will be red.
#' maptt(tt.events,
#'  multipleValuesByTime = c('infectedby', 'probabilities'),
#'  color1 = 'yellow',
#'  color2 = 'red',
#'  nbColors = 10,
#'  minWeight = 0,
#'  maxWeight = 30,
#'  weight1 = 20,
#'  weight2 = 25
#')
#'
#' @import htmlwidgets
#'
#' @export
maptt <- function(data,
                  multipleValuesByTime = c(),
                  circleRadius = 6,
                  defaultNodeColor = 'steelblue',
                  nodeColorByState = list(),
                  moveEdgeColor = 'steelblue',
                  color1 = 'green',
                  color2 = 'red',
                  nbColors = 10,
                  minWeight = 0,
                  maxWeight = 1,
                  weight1 = 0,
                  weight2 = 1,
                  autoFocus = TRUE,
                  keepOldFeatures = TRUE,
                  optionsControl = TRUE,
                  gradientControl = TRUE,
                  legend = TRUE,
                  width = NULL, height = NULL, elementId = NULL) {

  if (is.character(data) & length(data) == 1) {
    geoJson = data
  } else {
    geoJson = df2geojson(data, multipleValuesByTime)

    # Set default status color
    if (is.null(nodeColorByState) || length(nodeColorByState) == 0) {
      uniqueStatus <- unique(data$status[!is.na(data$status) & data$status != ''])
      nodeColorByState = createRainbowColors(uniqueStatus)
    }
  }

  args <- list(
    geoJson = geoJson,
    circleRadius = circleRadius,
    defaultNodeColor = defaultNodeColor,
    nodeColorByState = nodeColorByState,
    moveEdgeColor = moveEdgeColor,
    color1 = color1,
    color2 = color2,
    nbColors = nbColors,
    minWeight = minWeight,
    maxWeight = maxWeight,
    weight1 = weight1,
    weight2 = weight2,
    autoFocus = autoFocus,
    keepOldFeatures = keepOldFeatures,
    optionsControl = optionsControl,
    gradientControl = gradientControl,
    legend = legend
  )

  sizingPolicy <- htmlwidgets::sizingPolicy(
    viewer.padding = 0,
    viewer.fill = TRUE,
    browser.padding = 0,
    browser.fill = TRUE
  )

  htmlwidgets::createWidget("maptt", args, width = width, height = height, package = "SMITIDvisu", elementId = elementId, sizingPolicy = sizingPolicy)
}

#' mapttOutput
# Shiny bindings for maptt
#
# Output and render functions for using SMITIDvisu within Shiny
# applications and interactive Rmd documents.
#
# @param outputId output variable to read from
# @param width,height Must be a valid CSS unit (like \code{'100\%'},
#   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#   string and have \code{'px'} appended.
#'
#' @name SMITIDvisu-shiny
#'
#' @export
# Methods to make the widget work on Shiny
mapttOutput <- function(outputId, width = "100%", height = "400px") {
  shinyWidgetOutput(outputId, "maptt", width, height, package = "SMITIDvisu")
}


#' renderMaptt
# @param expr An expression that generates a SMITIDvisu
# @param env The environment in which to evaluate \code{expr}.
# @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#   is useful if you want to save an expression in a variable.
#' @rdname SMITIDvisu-shiny
#' @export
renderMaptt <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, mapttOutput, env, quoted = TRUE)
}
