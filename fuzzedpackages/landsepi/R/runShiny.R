#' runShinyApp
#' @description Launches landsepi shiny application into browser
#' @details R packages needed to run the shiny app :
#'  install.packages(c("shiny","DT", "shinyjs", "gridExtra", "png", "grid", "future", "promises", "tools"))
#' @importFrom utils installed.packages
#' @export
runShinyApp <- function() {
  appDir <- system.file("shiny-landsepi", "", package = "landsepi")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `landsepi`.", call. = FALSE)
  }
  
  needed_packages <- c("shiny","DT", "shinyjs", "gridExtra", "png", "grid", "future", "promises", "tools")
  if( sum(needed_packages %in% utils::installed.packages()[,1] == FALSE) != 0) {
    stop('Install packages : install.packages(c("shiny","DT", "shinyjs", "gridExtra", "png", "grid", "future", "promises", "tools"))')
  } 

  shiny::runApp(appDir, launch.browser = TRUE, display.mode = "normal")
}
