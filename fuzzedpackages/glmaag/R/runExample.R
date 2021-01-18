##' @title Shiny app
##' @description Run the shiny app for web interative using. You need to load data.table, dplyr, ggplot2, plotROC, and survminer beforehand. If you want to do parallel computing, you also need to register cores.
##' @param whetherrun whether to run shiny app, default to be TRUE
##' @importFrom shiny runApp
##' @examples runtheExample(FALSE)
##' @export
runtheExample <- function(whetherrun) {
  if (whetherrun) {
    appDir <- system.file("shiny_examples", "myapp.R", package = "glmaag")
    if (appDir == "") {
      stop("Could not find example directory. Try re-installing `glmaag`.", call. = FALSE)
    }
    runApp(appDir, display.mode = "normal")
  }
}