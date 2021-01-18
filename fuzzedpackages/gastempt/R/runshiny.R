#' Run shiny app demonstrating fit strategies with simulated data
#'
#' @return  Not used, starts shiny app
#' @export
#'
run_shiny = function() {
  appDir = system.file("shiny", package = "gastempt")
  if (appDir == "") {
    stop("Could not Shiny app in gastempt", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}