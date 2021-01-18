#Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=4113)' &

pkgs <- list()
pkgs[[ "shiny" ]] <- "1.0.0"
pkgs[[ "shinydashboard" ]] <- "0.5.3"
pkgs[[ "shinyjs" ]] <- "0.9"
pkgs[[ "ggplot2" ]] <- "2.2.1"
pkgs[[ "genepop" ]] <- "0.2.0"
pkgs[[ "yaml" ]] <- "2.1.14"
pkgs[[ "parallel" ]] <- "3.2.3"

for (package in c("shiny", "shinydashboard", "shinyjs", "ggplot2", "genepop", "yaml", "parallel")) {
  if (!require(package, character.only = T, quietly = T, warn.conflicts = FALSE)) {
    stop(paste("Package", package, "is not installed."))
  } else {
    vp <- packageVersion(package)
    if((pkgs[[package]] > vp)) {
      warning(paste("Version :" , vp, "of package", package, "is installed, genepop-shiny requires at least version:", pkgs[[package]]))
    }
  }
}

config <<- yaml.load_file("config.yml")$genepop
jscode <- "shinyjs.closeWindow = function() {open(location, '_self').close(); window.location.replace(document.URL); }"

set_restriction(TRUE)

source("./R/helper_functions.R")
source("./R/menugauche.R", local = T)
source("./pages/pages_def_home.R", local = T)
source("./pages/pages_def_load.R", local = T)
source("./pages/pages_def_opts1.R", local = T)
source("./pages/pages_def_opts2.R", local = T)
source("./pages/pages_def_opts3.R", local = T)
source("./pages/pages_def_opts4.R", local = T)
source("./pages/pages_def_opts5.R", local = T)
source("./pages/pages_def_opts6.R", local = T)
source("./pages/pages_def_opts7.R", local = T)
source("./pages/pages_def_opts8.R", local = T)
source("./pages/pages_def_documentation.R", local = T)

MODE_DEBUG <<- FALSE

style <- tags$style(HTML(readLines("www/added_styles.css")) )
UI <- dashboardPage(
  skin = config$skin,
  dashboardHeader(title=config$heading, titleWidth = 440, tags$li(hidden(div(id = "spinner", class = "plot-container", tags$img(src = config$spinner, id = "loading-spinner"))), class = "dropdown")),
  dashboardSidebar(width = 440, MenuGauche ),
  dashboardBody(
    shinyjs::useShinyjs(),
    extendShinyjs(text = jscode, functions = c("closeWindow")),
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.min.readable.css")) ,
    tags$head(style),
    tabItems(
      tabItem(tabName = "Home",         tabHome),
      tabItem(tabName = "Load",         tabLoad),
      tabItem(tabName = "Docu",         tabDocu),
      tabItem(tabName = "Opt11",        tabOpt11),
      tabItem(tabName = "Opt12",        tabOpt12),
      tabItem(tabName = "Opt13",        tabOpt13),
      tabItem(tabName = "Opt14",        tabOpt14),
      tabItem(tabName = "Opt15",        tabOpt15),
      tabItem(tabName = "Opt21",        tabOpt21),
      tabItem(tabName = "Opt22",        tabOpt22),
      tabItem(tabName = "Opt31",        tabOpt31),
      tabItem(tabName = "Opt32",        tabOpt32),
      tabItem(tabName = "Opt33",        tabOpt33),
      tabItem(tabName = "Opt34",        tabOpt34),
      tabItem(tabName = "Opt41",        tabOpt41),
      tabItem(tabName = "Opt51",        tabOpt51),
      tabItem(tabName = "Opt52",        tabOpt52),
      tabItem(tabName = "Opt53",        tabOpt53),
      tabItem(tabName = "Opt61",        tabOpt61),
      tabItem(tabName = "Opt62",        tabOpt62),
      tabItem(tabName = "Opt63",        tabOpt63),
      tabItem(tabName = "Opt64",        tabOpt64),
      tabItem(tabName = "Opt65",        tabOpt65),
      tabItem(tabName = "Opt66",        tabOpt66),
      tabItem(tabName = "Opt71",        tabOpt71),
      tabItem(tabName = "Opt72",        tabOpt72),
      tabItem(tabName = "Opt73",        tabOpt73),
      tabItem(tabName = "Opt74",        tabOpt74),
      tabItem(tabName = "Opt81",        tabOpt81),
      tabItem(tabName = "Opt82",        tabOpt82),
      tabItem(tabName = "Opt83",        tabOpt83),
      tabItem(tabName = "Opt84",        tabOpt84),
      tabItem(tabName = "Opt85",        tabOpt85),
      tabItem(tabName = "Opt86",        tabOpt86)
    )
  )
) #Interface


server <- function(input, output) {
  source("./opts/optLoad.R", local = TRUE)
  source("./opts/opt11.R", local = TRUE)
  source("./opts/opt12.R", local = TRUE)
  source("./opts/opt13.R", local = TRUE)
  source("./opts/opt14.R", local = TRUE)
  source("./opts/opt15.R", local = TRUE)
  source("./opts/opt21.R", local = TRUE)
  source("./opts/opt22.R", local = TRUE)
  source("./opts/opt31.R", local = TRUE)
  source("./opts/opt32.R", local = TRUE)
  source("./opts/opt33.R", local = TRUE)
  source("./opts/opt34.R", local = TRUE)
  source("./opts/opt41.R", local = TRUE)
  source("./opts/opt51.R", local = TRUE)
  source("./opts/opt52.R", local = TRUE)
  source("./opts/opt53.R", local = TRUE)
  source("./opts/opt61.R", local = TRUE)
  source("./opts/opt62.R", local = TRUE)
  source("./opts/opt63.R", local = TRUE)
  source("./opts/opt64.R", local = TRUE)
  source("./opts/opt65.R", local = TRUE)
  source("./opts/opt66.R", local = TRUE)
  source("./opts/opt71.R", local = TRUE)
  source("./opts/opt72.R", local = TRUE)
  source("./opts/opt73.R", local = TRUE)
  source("./opts/opt74.R", local = TRUE)
  source("./opts/opt81.R", local = TRUE)
  source("./opts/opt82.R", local = TRUE)
  source("./opts/opt83.R", local = TRUE)
  source("./opts/opt84.R", local = TRUE)
  source("./opts/opt85.R", local = TRUE)
  source("./opts/opt86.R", local = TRUE)

  observeEvent(input$close, {
    js$closeWindow()
    stopApp()
  })

  observeEvent(input$interrupt, {
    print("Interrupt")
  })

  DisableButtonGenepop()
}

shinyApp(ui = UI, server = server)
