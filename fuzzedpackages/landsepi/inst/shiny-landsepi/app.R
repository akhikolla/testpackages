library(shiny)
library(DT)
library(shinyjs)
#library(shinycssloaders)
library(landsepi)
data(package = "landsepi")

source("server.R")
source("ui.R")

if (interactive()) {
  shiny::shinyApp(ui = ui, server = server)
}
