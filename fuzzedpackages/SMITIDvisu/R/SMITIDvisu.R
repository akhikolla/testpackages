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

#' @encoding UTF-8
#' @title Visualize Data for Host and Viral Population from SMITIDstruct using HTMLwidgets
#'
#' @description Visualisation tools for SMITIDstruct package.
#' Allow to visualize host timeline, transmission tree, index diversities and variant graph using HTMLwidgets.
#' It mainly using D3JS, noUiSlider and FileSaver javascript libraries.
#'
#' @author Jean-Francois Rey \email{jean-francois.rey@@inrae.fr}
#' @author Julien Boge \email{julien.boge.u@@gmail.com}
#' @name SMITIDvisu-package
#' @aliases SMITIDvisu
#' @useDynLib SMITIDvisu, .registration=TRUE
#' @docType package
#' @details \tabular{ll}{
#'          Package: \tab SMITIDvisu\cr
#'          Type: \tab Package\cr
#'          Version: \tab 0.0.8\cr
#'          Date: \tab 2020-11-04\cr
#'          License: \tab GPL (>=3)\cr
#'          }
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @import yaml
#' @importFrom utils browseURL
#' @examples
#' \dontrun{
#'  library(SMITIDvisu)
#'  demo.SMITIDvisu.run()
#'}
NULL


#' Shiny bindings for visualisation widgets
#'
#' Output and render functions for using visualisation widgets within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{"100\%"},
#'   \code{"400px"}, \code{"auto"}) or a number, which will be coerced to a
#'   string and have \code{"px"} appended.
#' @param expr An expression that generates a networkD3 graph
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @importFrom htmlwidgets shinyWidgetOutput
#' @importFrom htmlwidgets shinyRenderWidget
#'
#' @name SMITIDvisu-shiny
NULL
