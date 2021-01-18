


tabOpt11 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt11-13.md"))))),
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
                             align="center",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$param,
                             sliderInput("Dememo11",label="Dememorization number", min=config$opts$o11$Dememo$min, max = config$opts$o11$Dememo$max, step=config$opts$o11$Dememo$step, value = config$opts$o11$Dememo$value),
                             sliderInput("Nbatches11", label="Number of batches", min=config$opts$o11$Nbatches$min, max = config$opts$o11$Nbatches$max, step=config$opts$o11$Nbatches$step, value = config$opts$o11$Nbatches$value),
                             sliderInput("Niters11", label="Number of iterations per batch", min=config$opts$o11$Niters$min, max = config$opts$o11$Niters$max, step=config$opts$o11$Niters$step, value = config$opts$o11$Niters$value),
                             actionButton("RunOpt11", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),

                           box(
                             title = "Dowload",
                             align="center",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$dowload,
                             style='height:235px;',
                             br(),
                             downloadButton('downloadOpt11All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,

                           box(
                             title="H1= Heterozygote deficiency: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             tabsetPanel(type = "tabs",
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Locus", verbatimTextOutput("Opt11outLoc")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Pops", verbatimTextOutput("Opt11outPop"))

                             )
                           ))

)

tabOpt12 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt11-13.md"))))),
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
                             align="center",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$param,
                             sliderInput("Dememo12",label="Dememorization number", min=config$opts$o12$Dememo$min, max = config$opts$o12$Dememo$max, step=config$opts$o12$Dememo$step, value = config$opts$o12$Dememo$value),
                             sliderInput("Nbatches12", label="Number of batches", min=config$opts$o12$Nbatches$min, max = config$opts$o12$Nbatches$max, step=config$opts$o12$Nbatches$step, value = config$opts$o12$Nbatches$value),
                             sliderInput("Niters12", label="Number of iterations per batch", min=config$opts$o12$Niters$min, max = config$opts$o12$Niters$max, step=config$opts$o12$Niters$step, value = config$opts$o12$Niters$value),
                             actionButton("RunOpt12", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:235px;',
                             br(),
                             downloadButton('downloadOpt12All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="H1= Heterozygote excess: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             tabsetPanel(type = "tabs",
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Locus", verbatimTextOutput("Opt12outLoc")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Pops", verbatimTextOutput("Opt12outPop"))
                             )
                           ))
)

tabOpt13 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt11-13.md"))))),
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
                             status = config$color$param,
                             align="center",
                             sliderInput("Dememo13",label="Dememorization number", min=config$opts$o13$Dememo$min, max = config$opts$o13$Dememo$max, step=config$opts$o13$Dememo$step, value = config$opts$o13$Dememo$value),
                             sliderInput("Nbatches13", label="Number of batches", min=config$opts$o13$Nbatches$min, max = config$opts$o13$Nbatches$max, step=config$opts$o13$Nbatches$step, value = config$opts$o13$Nbatches$value),
                             sliderInput("Niters13", label="Number of iterations per batch", min=config$opts$o13$Niters$min, max = config$opts$o13$Niters$max, step=config$opts$o13$Niters$step, value = config$opts$o13$Niters$value),
                             actionButton("RunOpt13", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             align="center",
                             status = config$color$dowload,
                             style='height:235px;',
                             br(),
                             downloadButton('downloadOpt13All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Probability test: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             tabsetPanel(type = "tabs",
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Locus", verbatimTextOutput("Opt13outLoc")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Pops", verbatimTextOutput("Opt13outPop")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By All Locus & All Pops", verbatimTextOutput("Opt13outLocPop"))
                             )
                           )
                     )
)

tabOpt14 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt14-15.md"))))),
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
                             sliderInput("Dememo14",label="Dememorization number", min=config$opts$o14$Dememo$min, max = config$opts$o14$Dememo$max, step=config$opts$o14$Dememo$step, value = config$opts$o14$Dememo$value),
                             sliderInput("Nbatches14", label="Number of batches", min=config$opts$o14$Nbatches$min, max = config$opts$o14$Nbatches$max, step=config$opts$o14$Nbatches$step, value = config$opts$o14$Nbatches$value),
                             sliderInput("Niters14", label="Number of iterations per batch",min=config$opts$o14$Niters$min, max = config$opts$o14$Niters$max, step=config$opts$o14$Niters$step, value = config$opts$o14$Niters$value),
                             actionButton("RunOpt14", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$dowload,
                             align="center",
                             style='height:235px;',
                             br(),
                             downloadButton('downloadOpt14All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Global Heterozygote deficiency: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             tabsetPanel(type = "tabs",
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Locus", verbatimTextOutput("Opt14outLoc")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Pops", verbatimTextOutput("Opt14outPop")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By All Locus & All Pops", verbatimTextOutput("Opt14outLocPop"))
                             )
                           ))
)
tabOpt15 = fluidRow(align="left",
                    box(
                      title="Help",
                      div(div(withMathJax(includeMarkdown(getFile("opt14-15.md"))))),
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
                             status = config$color$param,
                             align="center",
                             sliderInput("Dememo15",label="Dememorization number", min=config$opts$o15$Dememo$min, max = config$opts$o15$Dememo$max, step=config$opts$o15$Dememo$step, value = config$opts$o15$Dememo$value),
                             sliderInput("Nbatches15", label="Number of batches", min=config$opts$o15$Nbatches$min, max = config$opts$o15$Nbatches$max, step=config$opts$o15$Nbatches$step, value = config$opts$o15$Nbatches$value),
                             sliderInput("Niters15", label="Number of iterations per batch", min=config$opts$o15$Niters$min, max = config$opts$o15$Niters$max, step=config$opts$o15$Niters$step, value = config$opts$o15$Niters$value),
                             actionButton("RunOpt15", label = "Run", icon("paper-plane"), style=config$color$runButton)
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$dowload,
                             align="center",
                             style='height:235px;',
                             br(),
                             downloadButton('downloadOpt15All', label = "Results",  class=config$color$dowloadButton)
                           )),
                    column(width = 9,
                           box(
                             title="Global Heterozygote excess: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = config$color$result,
                             tabsetPanel(type = "tabs",
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Locus", verbatimTextOutput("Opt15outLoc")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Pops", verbatimTextOutput("Opt15outPop")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By All Locus & All Pops", verbatimTextOutput("Opt15outLocPop"))
                             )
                           ))
)
