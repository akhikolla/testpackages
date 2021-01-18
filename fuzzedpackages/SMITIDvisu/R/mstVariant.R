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


#' mstVariant
#' @description Draw Variants genotypes distances as a graph using Minimum Spanning Tree algorithm.
#' @param mat a distance matrix between sequence of variants (interger distance no floating values)
#' @param prop a data.frame for variants sequences proportions and count (see details)
#' @param node.prop list of variants with proportions and time (default NULL)
#' @param width numeric width for the area in pixels.
#' @param height numeric hieght for the area in pixels.
#' @param elementId the element ID where is draw
#' 
#' @details
#' \strong{mat} is a simple distance matrix with interger values, row and lines contain a unique identifier of each variant sequences.
#' \strong{prop} is a data.frame where each row is a variant sequence, it have to contain in columns factor "ID", "proportion" and "count".
#' "ID" is a unique identifier matching matrix value identifier, "proportion" is the proportions of the variant sequence
#' and "count" the number of variant sequence in a varions set.
#' \strong{node.prop} is a list with name that matching \strong{mat} identifier and \strong{prop} "ID". Each list element contains a subvector time (Julian or timestamp) and value (proportions).
#' That allow to draw variants proportions over time.
#'
#' @examples 
#' library(SMITIDvisu)
#' data(st)
#' mstV <- mstVariant(st.dist113_all,st.prop113_all, st.listTimeProp113)
#' \dontrun{ 
#' ## export as standalone html file
#' htmlwidgets::saveWidget(mstV, "mstVariant.html")
#' browseURL("mstVariant.html")
#' }
#'
#' @import htmlwidgets
#'
#' @export
mstVariant <- function(mat,
                       prop,
                       node.prop = NULL,
                       width = NULL,
                       height = NULL,
                       elementId = NULL) {
  
  
  # transform the distance matrix and varions proportions into a list
  graph.list = c(createMSTGraph(mat,prop), "node_prop"=list(unname(node.prop)))

  # forward options using x
  x = list(
    data = jsonlite::toJSON(graph.list),
    options = list()
  )
  
  # create widget
  htmlwidgets::createWidget(
    name = 'mstVariant',
    x,
    width = width,
    height = height,
    package = 'SMITIDvisu',
    elementId = elementId
  )
}

#' mstVariantOutput
# Shiny bindings for SMITIDvisu
#
# Output and render functions for using SMITIDvisu within Shiny
# applications and interactive Rmd documents.
#
# @param outputId output variable to read from
# @param width,height Must be a valid CSS unit (like \code{'100\%'},
#   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#   string and have \code{'px'} appended.
# @param expr An expression that generates a SMITIDvisu
# @param env The environment in which to evaluate \code{expr}.
# @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#   is useful if you want to save an expression in a variable.
#
#' @name SMITIDvisu-shiny
#'
#' @export
mstVariantOutput <- function(outputId, width = '100%', height = '600px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'mstVariant', width, height, package = 'SMITIDvisu')
}

#' rendermstVariant
#
# @param expr An expression that generates a SMITIDvisu
# @param env The environment in which to evaluate \code{expr}.
# @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#   is useful if you want to save an expression in a variable.
#'
#' @rdname SMITIDvisu-shiny
#' @export
rendermstVariant <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, mstVariantOutput, env, quoted = TRUE)
}
