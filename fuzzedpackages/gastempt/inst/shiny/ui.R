# Libraries are included to keep intellisense quiet
suppressPackageStartupMessages(library(shinyjs))
library(shinyAce)
library(shinyBS)

shinyUI(
  fluidPage(
    useShinyjs(),
    theme = "bootstrap.css",
    tags$link(rel = "stylesheet", type = "text/css", href = "gastempt.css"),
    titlePanel("Fitting gastric emptying curves"),
    sidebarLayout(
      sidebarPanel(
        h3("Analyze data"),
        selectInput("method_a", "Method", choices =
                      c("nlme population fit" = "nlme",
                       "Bayesian Stan fit" = "stan"),
                    selected = "nlme"),
        selectInput("fit_model", "Curves analyzed as ", choices =
                      c("linexp (with overshoot)" = "linexp",
                        "powexp (without overshoot)" = "powexp"),
                    selected = "linexp"),
        conditionalPanel(
          condition = "input.method_a == 'nlme'",
          selectInput("variant", "Variant", choices =
                    c("1: v0 + tempt + kappa" = 1,
                      "2: pdDiag(v0 + tempt + kappa)" = 2,
                      "3: v0 + tempt" = 3))
        ),
        conditionalPanel(
          condition = "input.method_a == 'stan'",
          selectInput("cov_model", "Model name", choices =
                      c("No covariance" = "nocov",
                        "With covariance" = "withcov")),
          conditionalPanel(
            condition = "input.manual",
            fluidRow(splitLayout(
              numericInput("lkj", "LKJ", 2, 1, 50, 1, width = "80px"),
              numericInput("student_df", "Student df", 5, 2, 30, 1, width = "80px")
            ))
          )
        ), # conditionalPanel stan
        checkboxInput("manual", "Expert Mode"),
        hr(),
        h3("Simulated data"),
        h4("Preset"),
        selectInput("preset",NULL,
                    choices = setNames(presets$id, presets$label)),
        conditionalPanel(
          condition = "input.manual",
          h4("Manual settings"),
          fluidRow(splitLayout(
            numericInput("n_records", "# records", 12, 5, 30, 5, width = "80px"),
            numericInput("seed", "Seed", 33, width = "80px")
          )),
          selectInput("model_s","Data generated as ", choices =
                        c("linexp (with overshoot)" = "linexp",
                          "powexp (without overshoot)" = "powexp")),
          fluidRow(splitLayout(
            h5(HTML("<br><br>v0")),
            numericInput("v0_mean", "Mean", 400,  50, 800, 50),
            numericInput("v0_std_perc", "Std (%)", 20, 5, 30, 5),
            cellWidths = c("25%","30%","30%"))),
          fluidRow(splitLayout(
            h5("tempt"),
            numericInput("tempt_mean", NULL, 60, 20 , 200, 5),
            numericInput("tempt_std_perc", NULL, 25, 5, 30, 5),
            cellWidths = c("25%","30%","30%"))),
          fluidRow(splitLayout(
            h5(HTML("beta/<br>kappa")),
            numericInput("kappa_beta_mean", NULL, 1, 0.1, 4, 0.1),
            numericInput("kappa_beta_std_perc", NULL, 10, 5, 30, 5),
            cellWidths = c("25%","30%","30%"))),

          div(class = 'withSeparator',
            selectInput("student_t_df","Noise type/outliers",
                        choices = c("Gauss" = 0,
                                    "weak (Student-7)" = 7,
                                    "some (Student-5)" = 5,
                                    "heavy (Student-3)" = 3))),
          fluidRow(splitLayout(
            numericInput("noise_perc", "Noise %", 20, 5, 50, 5, width = "60px"),
            numericInput("missing", "Missing %", 0, 0, 50, 5, width = "60px")
          )),
          hr()
        ), # conditionalPanel manual,
        helpText(HTML('<a href = "https://github.com/dmenne/gastempt" target="_blank">github source code</a>')),
        helpText(HTML('<a href = "https://github.com/dmenne/gastempt/issues" target="_blank">Report issues</a>')),
        checkboxInput("show_pop", "Show popover help", value = TRUE),
        # The following should not be moved to the server, because one
        # popover must remain on the page to load dependencies
        bsPopover("show_pop",  "Enable/disable all popups","","right"),
        width = 3),
      mainPanel(
        tabsetPanel(
          tabPanel("Data",
            aceEditor("data", "", mode = "plain_text"),
            actionButton("clearButton","Clear", icon = icon("eraser")),
            tags$script(type = "text/javascript",
              HTML("ace.edit('data').setOptions({tabSize:12,showInvisibles:true,useSoftTabs:false});")),
            div(id = "plot-container",
                tags$img(src = "spinner.gif",
                         id = "loading-spinner"),
                plotOutput("fit_plot")
            ),
            hr(),
            downloadButton("download_coef", "Download"),
            DT::dataTableOutput("table")
          ), # tabPanel
          tabPanel("Check",
              plotOutput("residual_plot"),
              plotOutput("trace_v")
            ) # check tabPanel
        ) # tabsetPanel
      ) # mainPanel
    ) # sidebarLayout
  ) # fluidPanel
) # shinyUI

