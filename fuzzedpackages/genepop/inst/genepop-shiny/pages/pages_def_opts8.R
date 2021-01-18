tabOpt81 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt81.md"))))),
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
                             radioButtons("methodAllele81", label = h3("Null Allele Method"),
                                          choices = list("Default" = "Default",
                                                         "ApparentNulls" = "ApparentNulls", "B96" = "B96"), selected = "Default"),
                             sliderInput("coverage81",label="CIcoverage", min=0, max = 1, step=0.01, value = 0.95),
                             actionButton("RunOpt81", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt81All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Allele and genotype frequencies",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt81out")
                           ))
)

tabOpt82 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt82.md"))))),
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
                             actionButton("RunOpt82", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt82All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Diploidisation of haploid data",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt82out")
                           ))
)

tabOpt83 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt83.md"))))),
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
                             actionButton("RunOpt83", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt83All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Relabeling alleles",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt83out")
                           ))
)

tabOpt84 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt84-85.md"))))),
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
                             actionButton("RunOpt84", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt84All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Conversion to individual data with population names",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt84out")
                           ))
)

tabOpt85 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt84-85.md"))))),
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
                             actionButton("RunOpt85", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt85All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Conversion to individual data with individual names",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt85out")
                           ))
)

tabOpt86 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt86.md"))))),
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
                             actionButton("RunOpt86", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt86All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Random sampling of haploid genotypes from diploid ones",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt86out")
                           ))
)
