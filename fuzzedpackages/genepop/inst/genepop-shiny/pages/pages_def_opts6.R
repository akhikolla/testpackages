tabOpt61 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt61-64.md"))))),
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
                             radioButtons("ploidy61", label = h3("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             actionButton("RunOpt61", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt61All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Allele identity (F-statistics) for all populations",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt61out")
                           ))
)

tabOpt62 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt61-64.md"))))),
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
                             radioButtons("ploidy62", label = h3("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             actionButton("RunOpt62", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt62All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Allele identity (F-statistics) for all population pairs",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt62out")
                           ))
)

tabOpt63 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt61-64.md"))))),
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
                             radioButtons("ploidy63", label = h3("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             actionButton("RunOpt63", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt63All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Allele size (Rho-statistics) for all populations",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt63out")
                           ))
)

tabOpt64 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt61-64.md"))))),
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
                             radioButtons("ploidy64", label = h3("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             actionButton("RunOpt64", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt64All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Allele size (Rho-statistics) for all population pairs",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt64out")
                           ))
)

tabOpt65 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt65.md"))))),
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
                             radioButtons("ploidy65", label = h4("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             radioButtons("isolationStatistic65", label = h4("Isolation Statistic"),
                                          choices = list("a" = "a", "e" = "e"), selected = "e"),
                             radioButtons("geographicScale65", label = h4("Geographic Scale"),
                                          choices = list("Log" = "2D", "Linear" = "1D"), selected = "2D"),
                             sliderInput("coverage65",label="CIcoverage", min=0, max = 1, step=0.01, value = 0.95),
                             numericInput("minDistance65", "Minimal geographic distance", 0.0001, min = 0, max = 1000000000),
                             numericInput("maxDistance65", "Maximal geographic distance", 1000000000, min = 0, max = 1000000000),
                             textInput("testPoint65", "Test Point", "NA"),
                             numericInput("mantelSeed65", "Mantel Seed", 87654321, min = 0, max = 1000000000),
                             radioButtons("mantelRankTest65", label = h4("Mantel Rank Test"),
                                          choices = list("True" = "True", "False" = "False"), selected = "False"),
                             numericInput("mantelPermutations65", "Mantel Permutations", 1000, min = 0, max = 100000),
                             actionButton("RunOpt65", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt65All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Isolation by distance between individuals",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt65out")
                           ))
)

tabOpt66 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt66.md"))))),
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
                             radioButtons("ploidy66", label = h4("Ploidy"),
                                          choices = list("Haploid" = "Haploid", "Diploid" = "Diploid"), selected = "Diploid"),
                             radioButtons("isolationStatistic66", label = h4("Isolation Statistic"),
                                          choices = list("Fst/(1-Fst)" = "F/(1-F)", "Common denominator" = "SingleGeneDiv"), selected = "F/(1-F)"),
                             radioButtons("geographicScale66", label = h4("Geographic Scale"),
                                          choices = list("Log" = "2D", "Linear" = "1D"), selected = "2D"),
                             sliderInput("coverage66",label="CIcoverage", min=0, max = 1, step=0.01, value = 0.95),
                             numericInput("minDistance66", "Minimal geographic distance", 0.0001, min = 0, max = 1000000000),
                             numericInput("maxDistance66", "Maximal geographic distance", 1000000000, min = 0, max = 1000000000),
                             textInput("testPoint66", "Test Point", "NA"),
                             numericInput("mantelSeed66", "Mantel Seed", 87654321, min = 0, max = 1000000000),
                             radioButtons("mantelRankTest66", label = h4("Mantel Rank Test"),
                                          choices = list("True" = "True", "False" = "False"), selected = "False"),
                             numericInput("mantelPermutations66", "Mantel Permutations", 1000, min = 0, max = 100000),
                             actionButton("RunOpt66", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:100px;',
                             br(),
                             downloadButton('downloadOpt66All', label = "Results",  class=config$color$dowloadButton)
                           )
                    ),
                    column(width = 9,
                           box(
                             title="Isolation by distance between groups",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             style='height:600px; overflow-y: scroll',
                             verbatimTextOutput("Opt66out")
                           ))
)
