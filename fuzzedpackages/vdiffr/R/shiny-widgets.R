#' Shiny bindings for graphical comparison widgets
#'
#' Output and render functions for using the Toggle, Slide and Diff
#' widgets within Shiny applications and interactive Rmd documents.
#' Used in [manage_cases()].
#'
#' @param outputId Output variable to read from.
#' @param width,height Must be a valid CSS unit (like `"100\%"`,
#'   `"400px"`, `"auto"`) or a number, which will be coerced to a
#'   string and have `"px"` appended.
#' @param expr An expression that generates a comparison widget.
#' @param env The environment in which to evaluate `expr`.
#' @param quoted Is `expr` a quoted expression (with `quote()`)? This
#'   is useful if you want to save an expression in a variable.
#' @name shinybindings
NULL

#' @rdname shinybindings
#' @export
toggleOutput <- function(outputId, width = "100%", height = "400px") {
  htmlwidgets::shinyWidgetOutput(outputId, "vdiffr-toggle",
    width, height, package = "vdiffr")
}

#' @rdname shinybindings
#' @export
renderToggle <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) {
    expr <- substitute(expr)
  }
  htmlwidgets::shinyRenderWidget(expr, toggleOutput, env, quoted = TRUE)
}

#' @rdname shinybindings
#' @export
slideOutput <- function(outputId, width = "100%", height = "400px") {
  htmlwidgets::shinyWidgetOutput(outputId, "vdiffr-slide",
    width, height, package = "vdiffr")
}

#' @rdname shinybindings
#' @export
slideTransition <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) {
    expr <- substitute(expr)
  }
  htmlwidgets::shinyRenderWidget(expr, slideOutput, env, quoted = TRUE)
}

#' @rdname shinybindings
#' @export
diffOutput <- function(outputId, width = "100%", height = "400px") {
  htmlwidgets::shinyWidgetOutput(outputId, "vdiffr-diff",
    width, height, package = "vdiffr")
}

#' @rdname shinybindings
#' @export
diffTransition <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) {
    expr <- substitute(expr)
  }
  htmlwidgets::shinyRenderWidget(expr, diffOutput, env, quoted = TRUE)
}
