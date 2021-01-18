#' timeLine
#' @description Draw a host time line.
#' Time use timestamp or Date in ISO format.
#' @param data a data.frame that reprensent hosts status in time with ID, status and time in columns
#' @param title a tttle as character
#' @param color list of color for timeline elements
#' @param width numeric width for the area in pixels.
#' @param height numeric hieght for the area in pixels.
#' @param elementId the element ID where is draw
#'
#' @examples 
#' library(SMITIDvisu)
#' data(hostline)
#' tl <- timeLine(hostline,
#'                title="Example host 113",
#'                color=list("infected"="red","offspring"="green",
#'                              "alive"="blue","inf"="orange",
#'                              "dead"="black","Obs"="purple"))
#' \dontrun{
#' ## export as standalone html file
#' htmlwidgets::saveWidget(tl, "timeline.html")
#' browseURL("timeline.html")
#' }
#'
#' @import htmlwidgets
#' @import magrittr
#'
#' @export
timeLine <- function(data,
                     title,
                     color = NULL,
                     width = NULL,
                     height = NULL,
                     elementId = NULL) {

  if( is.null(color)) {
    color = createRainbowColors(unique(data$label))
  }
  
  # forward options using x
  x = list(
    data = jsonlite::toJSON(createTimeLine(data, title)),
    #data = createTimeLine(data, title),
    options = list(color=color)
  )
  
  # create widget
  htmlwidgets::createWidget(
    name = 'timeLine',
    x,
    width = width,
    height = height,
    package = 'SMITIDvisu',
    elementId = elementId
  )
}

#' timeLineOuput 
# Shiny bindings for timeline
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
timeLineOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'timeLine', width, height, package = 'SMITIDvisu')
}

#' renderTimeLine
#
# @param expr An expression that generates a SMITIDvisu
# @param env The environment in which to evaluate \code{expr}.
# @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#   is useful if you want to save an expression in a variable.
#'
#' @rdname SMITIDvisu-shiny
#' @export
renderTimeLine <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, timeLineOutput, env, quoted = TRUE)
}
