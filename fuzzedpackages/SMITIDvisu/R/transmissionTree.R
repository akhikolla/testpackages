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


#' transmissionTree
#' @description Draw a transmission tree over the time.
#' Time use timestamp or Date in ISO format ("%Y-%m-%dT%H:%M:%S").
#' @param nodes a data.frame that reprensent hosts status in time with ID, status and time in columns
#' @param edges a data.frame that reprensent tramsmission link between hosts (pathogens) with ID, source, target and time in columns
#' @param nodes.color a list of color for nodes status "status"="color"
#' @param width numeric width for the area in pixels.
#' @param height numeric hieght for the area in pixels.
#' @param elementId the element ID where is draw
#'
#' @examples 
#' library(SMITIDvisu)
#' data(transmissiontree)
#' tt <- transmissionTree(tt.nodes,tt.edges, nodes.color = list("default"="black","Inf"="red"))
#' \dontrun{
#' ## export as standalone html file
#' htmlwidgets::saveWidget(tt, "transTree.html")
#' browseURL("transTree.html")
#' }
#'
#' @import htmlwidgets
#'
#' @export
transmissionTree <- function(nodes,
                             edges, 
                             nodes.color = NULL,
                             width = NULL,
                             height = NULL,
                             elementId = NULL) {

  
  
  graph.list = createTimeGraph(nodes,edges)
  
  if( is.null(nodes.color)) {
    nodes.color = createRainbowColors(unique(nodes$status))
  }
  
  # forward options using x
  x = list(
    data = jsonlite::toJSON(graph.list),
    options = list(nodes_color=nodes.color)
  )
  
  # create widget
  htmlwidgets::createWidget(
    name = 'transmissionTree',
    x,
    width = width,
    height = height,
    package = 'SMITIDvisu',
    elementId = elementId
  )
}

#' transmissionTreeOutput
# Shiny bindings for transmissionTree
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
transmissionTreeOutput <- function(outputId, width = '100%', height = "500px"){
  htmlwidgets::shinyWidgetOutput(outputId, 'transmissionTree', width, height, package = 'SMITIDvisu')
}

#' renderTransmissionTree
# @param expr An expression that generates a SMITIDvisu
# @param env The environment in which to evaluate \code{expr}.
# @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#   is useful if you want to save an expression in a variable.
#' @rdname SMITIDvisu-shiny
#' @export
renderTransmissionTree <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, transmissionTreeOutput, env, quoted = TRUE)
}
