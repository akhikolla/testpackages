library(shiny)
library(DT)
library(shinyjs)

# UI
######################################################################################
######################################################################################
landscapeTab <- {
  shiny::tabPanel(
    "Landscape",
    shiny::br(),
    shiny::fluidRow(
      column(
        width = 6,
        shiny::selectInput(
          inputId = "landscape",
          label = "Landscape structure (field boundaries)",
          choices = list(
            "Landscape 1" = 1,
            "Landscape 2" = 2,
            "Landscape 3" = 3,
            "Landscape 4" = 4,
            "Landscape 5" = 5
          ),
          selected = 1,
        )
      ),
      column(
        width = 6,
        shiny::selectInput(
          inputId = "aggregLevel",
          label = "Spatial aggregation of croptypes",
          choices = list(
            "Highly fragmented" = "low",
            "Balanced" = "medium",
            "Highly aggregated" = "high"
          ),
          selected = "low",
        )
      )
    ),
    hr(),
    shiny::fluidRow(
      tags$div(lang="en",
               column(
                 width = 3,
                 PercentageInput(
                   inputId = "prop0",
                   label = "C0 proportion",
                   value =  0.33
                 )
               ),
               column(
                 width = 3,
                 PercentageInput(
                   inputId = "prop1",
                   label = "C1 proportion",
                   value =  0.33
                 )
               ),
               column(
                 width = 3,
                 PercentageInput(
                   inputId = "prop2",
                   label = "C2 proportion",
                   value =  0.34
                 )
               )
      ),
      column(
        width = 3,
        IntegerInput(
          inputId = "rotationPeriod",
          label = "Rotation period (years)",
          value = 1,
          max = 50
        )
      )
    ),
    hr(),
    shiny::fluidRow(
      column(
        width = 4,
        IntegerInput(
          inputId = "nYear",
          label = "Simulation duration (years)",
          value = 30,
          max = 50
        )
      ),
      column(
        width = 4,
        IntegerInput(
          inputId = "nTSpY",
          label = "Time steps per year (days)",
          value = 120,
          max = 365
        )
      ),
      column(
        width = 4,
        IntegerInput(
          inputId = "seed",
          label = "Seed (for random number generator)",
          value = 12345,
          max = 99999
        )
      )
    )
  )
}
######################################################################################
cultivarTab <- {
  shiny::tabPanel(
    "Croptypes and Cultivars",
    DT::DTOutput(outputId = "croptypes"),
    hr(),
    DT::DTOutput(outputId = "cultivars")
  )
}
######################################################################################
inputUi <- {
  shiny::sidebarPanel(
    shiny::h3("Input"),
    shiny::div(
      shiny::selectInput(
        inputId = "demo",
        label = "Default Strategies",
        choices = list(
          "Mosaic" = "MO",
          "Mixture" = "MI",
          "Rotation" = "RO",
          "Pyramiding" = "PY"
        ),
        width = "25%"
      ),
      align = "center"
    ),
    shiny::tabsetPanel(landscapeTab, cultivarTab),
    width = 6,
    align = "center",
    id = "inputpanel"
  )
}
######################################################################################
outputUi <- {
  shiny::mainPanel(
    shiny::h3("Output"),
    shiny::plotOutput(outputId = "landscapeimg", dblclick = "plot_landscapeimg"),
    shiny::uiOutput(outputId = "video"),
    width = 6,
    align = "center"
  )
}
######################################################################################
ui <- {
  shiny::fluidPage(
    tags$html(lang="en"),
    title = "Landsepi Demo",
    shinyjs::useShinyjs(),
    tags$head(
      tags$link(href = "style.css", rel = "stylesheet"),
      tags$script(src = "script.js"),
      tags$meta(name = "description", content = "A stochastic, spatially-explicit, demo-genetic model 
                simulating the spread and evolution of a plant pathogen in a heterogeneous landscape 
                to assess resistance deployment strategies. It is based on a spatial geometry for describing 
                the landscape and allocation of different cultivars, a dispersal kernel for the 
                dissemination of the pathogen, and a SEIR ('Susceptible-Exposed-Infectious-Removed’) 
                structure with a discrete time step. It provides a useful tool to assess the performance 
                of a wide range of deployment options with respect to their epidemiological, 
                evolutionary and economic outcomes."),
      tags$meta(name = "title", content = "Landsepi: Landscape Epidemiology and Evolution"),
      tags$meta(name = "author", content = "Jean-François Rey"),
      tags$meta(name = "keywords", content = "R, Shiny, INRA, landscape, epidemiology, strategies, deployment, hosts, pathogens, stochastic, spatio-temporal, demonstration")
    ),
    fluidRow(
    titlePanel("Landsepi : Landscape Epidemiology and Evolution"),
    actionButton("About", "About"),
    align="center"),
    shiny::br(),
    shiny::sidebarLayout(inputUi, outputUi),
    shiny::fluidRow(shiny::div(
      shiny::actionButton(inputId = "generateLandscape", label = "Generate the landscape"),
      align = "center"
    )),
    shiny::br(),
    shiny::fluidRow(shiny::div(
      shiny::actionButton(inputId = "runSimulation", label = "Run simulation"),
      shiny::actionButton(inputId = "stopSimulation", label = "Stop simulation"),
      align = "center"
    )),
    shiny::br(),
    shiny::fluidRow(shiny::div(
      shiny::downloadButton(outputId = "export", label = "Export simulation"),
      align = "center"
    )),
    shiny::br()
  )
}
