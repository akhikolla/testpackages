tabOpt51 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt51.md"))))),
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
                             actionButton("RunOpt51", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt51All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Allele and genotype frequencies per locus and per sample",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt51out")
                           ))
)

tabOpt52 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt52.md"))))),
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
                             radioButtons("ploidy52", label = h3("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             actionButton("RunOpt52", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt52All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Gene diversities and Fis using allele identity",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt52out")
                           ))
)

tabOpt53 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt53.md"))))),
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
                             radioButtons("ploidy53", label = h3("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             actionButton("RunOpt53", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt53All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Gene diversities and Fis using allele size",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt53out")
                           ))
)
