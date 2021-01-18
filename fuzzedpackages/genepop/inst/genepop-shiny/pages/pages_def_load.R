tabLoad = fluidPage(align="left",
box(
  title="Help",
  div(HTML("Choose a file from you local disc to upload on the server. </br>The content of your file will be summarized and loaded in a table")),
  width = 12,
  solidHeader = FALSE,
  status = config$color$help,
  style='height:150px; overflow-y: scroll',
  collapsible = TRUE,
  collapsed = TRUE
),
                    box(
                      title = "Input file",
                      width = 4,
                      solidHeader = TRUE,
                      status = config$color$input,
                      fileInput('file1', 'Upload Genepop File',
                                accept=NULL),
                      #plotOutput("distPops"),
                      numericInput("randomSeed", "Random Seed", 12345678, min = 1, max = 1000000000)
                    ),
                    box(
                      title = "File content",
                      width = 8,
                      solidHeader = TRUE,
                      status = config$color$result,
                      tabsetPanel(type = "tabs",
                                  tabPanel("Tab", div(style = 'overflow-x: scroll', DT::dataTableOutput('contents'))),
                                  tabPanel(style='height:600px; overflow-y: scroll; overflow-x: scroll', "Input file", verbatimTextOutput("InputFileData")))
                    )

)
