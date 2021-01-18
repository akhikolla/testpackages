tabOpt21 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt21.md"))))),
                      width = 12,
                      solidHeader = FALSE,
                      status = "warning",
                      style='height:250px; overflow-y: scroll',
                      collapsible = TRUE,
                      collapsed = TRUE
                    ),
                    column(width = 3,
                           box(
                             title = "Params",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$param,
                             sliderInput("Dememo21",label="Dememorization number", min=config$opts$o21$Dememo$min, max = config$opts$o21$Dememo$max, step=config$opts$o21$Dememo$step, value = config$opts$o21$Dememo$value),
                             sliderInput("Nbatches21", label="Number of batches", min=config$opts$o11$Nbatches$min, max = config$opts$o11$Nbatches$max, step=config$opts$o11$Nbatches$step, value = config$opts$o11$Nbatches$value),
                             sliderInput("Niters21", label="Number of iterations per batch", min=config$opts$o21$Niters$min, max = config$opts$o21$Niters$max, step=config$opts$o21$Niters$step, value = config$opts$o21$Niters$value),
                             actionButton("RunOpt21", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt21All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="LD for each pair of loci in each population: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt21out")
                           ))
)

tabOpt22 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt22.md"))))),
                      width = 12,
                      solidHeader = FALSE,
                      status = "warning",
                      style='height:250px; overflow-y: scroll',
                      collapsible = TRUE,
                      collapsed = TRUE
                    ),
                    column(width = 3,
                           box(
                             title = "Params",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$param,
                             actionButton("RunOpt22", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt22All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="genotypic contingency tablesgenotypic contingency tables: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt22out")
                           ))
)
