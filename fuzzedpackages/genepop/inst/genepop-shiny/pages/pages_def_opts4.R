tabOpt41 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt41.md"))))),
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
                             radioButtons("ploidy41", label = h3("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             actionButton("RunOpt41", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt41All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Nm estimates (private allele method): Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt41out")
                           ))
)
