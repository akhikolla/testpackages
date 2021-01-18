tabOpt31 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt31-32.md"))))),
                      width = 12,
                      solidHeader = FALSE,
                      status = config$color$help,
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
                             sliderInput("Dememo31",label="Dememorization number", min=config$opts$o31$Dememo$min, max = config$opts$o31$Dememo$max, step=config$opts$o31$Dememo$step, value = config$opts$o31$Dememo$value),
                             sliderInput("Nbatches31", label="Number of batches", min=config$opts$o31$Nbatches$min, max = config$opts$o31$Nbatches$max, step=config$opts$o31$Nbatches$step, value = config$opts$o31$Nbatches$value),
                             sliderInput("Niters31", label="Number of iterations per batch", min=config$opts$o31$Niters$min, max = config$opts$o31$Niters$max, step=config$opts$o31$Niters$step, value = config$opts$o31$Niters$value),
                             actionButton("RunOpt31", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt31All', label = "Results",  class="btn-primary")
                           )
                    ),column(width = 9,
                             box(
                               title="Genic differentiation for all populations : Results",
                               width = NULL,
                               solidHeader = TRUE,
                               status = config$color$result,
                               style='height:600px; overflow-y: scroll',
                               verbatimTextOutput("Opt31out")
                             ))
)

tabOpt32 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt31-32.md"))))),
                      width = 12,
                      solidHeader = FALSE,
                      status = config$color$help,
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
                             sliderInput("Dememo32",label="Dememorization number", min=config$opts$o32$Dememo$min, max = config$opts$o32$Dememo$max, step=config$opts$o32$Dememo$step, value = config$opts$o32$Dememo$value),
                             sliderInput("Nbatches32", label="Number of batches", min=config$opts$o32$Nbatches$min, max = config$opts$o32$Nbatches$max, step=config$opts$o32$Nbatches$step, value = config$opts$o32$Nbatches$value),
                             sliderInput("Niters32", label="Number of iterations per batch", min=config$opts$o32$Niters$min, max = config$opts$o32$Niters$max, step=config$opts$o32$Niters$step, value = config$opts$o32$Niters$value),
                             actionButton("RunOpt32", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt32All', label = "Results",  class="btn-primary")
                           )
                    ),column(width = 9,
                             box(
                               title="Genic differentiation for all pairs of populations : Results",
                               width = NULL,
                               solidHeader = TRUE,
                               status = config$color$result,
                               style='height:600px; overflow-y: scroll',
                               verbatimTextOutput("Opt32out")
                             ))
)

tabOpt33 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt33-34.md"))))),
                      width = 12,
                      solidHeader = FALSE,
                      status = config$color$help,
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
                             sliderInput("Dememo33",label="Dememorization number", min=config$opts$o33$Dememo$min, max = config$opts$o33$Dememo$max, step=config$opts$o33$Dememo$step, value = config$opts$o33$Dememo$value),
                             sliderInput("Nbatches33", label="Number of batches", min=config$opts$o33$Nbatches$min, max = config$opts$o33$Nbatches$max, step=config$opts$o33$Nbatches$step, value = config$opts$o33$Nbatches$value),
                             sliderInput("Niters33", label="Number of iterations per batch", min=config$opts$o33$Niters$min, max = config$opts$o33$Niters$max, step=config$opts$o33$Niters$step, value = config$opts$o33$Niters$value),
                             actionButton("RunOpt33", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt33All', label = "Results",  class="btn-primary")
                           )
                    ),column(width = 9,
                             box(
                               title="Genotypic differentiation for all populations : Results",
                               width = NULL,
                               solidHeader = TRUE,
                               status = config$color$result,
                               style='height:600px; overflow-y: scroll',
                               verbatimTextOutput("Opt33out")
                             ))
)

tabOpt34 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt33-34.md"))))),
                      width = 12,
                      solidHeader = FALSE,
                      status = config$color$help,
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
                             sliderInput("Dememo34",label="Dememorization number", min=config$opts$o34$Dememo$min, max = config$opts$o34$Dememo$max, step=config$opts$o34$Dememo$step, value = config$opts$o34$Dememo$value),
                             sliderInput("Nbatches34", label="Number of batches", min=config$opts$o34$Nbatches$min, max = config$opts$o34$Nbatches$max, step=config$opts$o34$Nbatches$step, value = config$opts$o34$Nbatches$value),
                             sliderInput("Niters34", label="Number of iterations per batch", min=config$opts$o34$Niters$min, max = config$opts$o34$Niters$max, step=config$opts$o34$Niters$step, value = config$opts$o34$Niters$value),
                             actionButton("RunOpt34", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt34All', label = "Results",  class="btn-primary")
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Genotypic differentiation for all pairs of populations : Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt34out")
                           ))
)
