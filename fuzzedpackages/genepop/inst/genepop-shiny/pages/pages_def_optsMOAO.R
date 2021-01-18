tabOptMOAO1 = fluidRow(align="left",
                    box(
                      title="Help",
                      tags$iframe(style="height:600px; width:100%",frameborder="0",scrolling="no", src="http://kimura.univ-montp2.fr/~rousset/Genepop.pdf#subsubsection.4.1.3"),
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
                             status = "primary",
                             fileInput('fileMOAO1', 'Upload ad-hoc HW file', accept=NULL),
                             sliderInput("DememoMOAO1",label="Dememorization number", min=100, max = 100000, step=10, value = 10000),
                             sliderInput("NbatchesMOAO1", label="Number of batches", min = 10, max = 10000, step=10, value = 100),
                             sliderInput("NitersMOAO1", label="Number of iterations per batch", min = 400, max = 100000, step=10, value = 5000),
                             actionButton("RunOptMOAO1", label = "Run", icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                           ),

                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             status = "primary",
                             style='height:235px;',
                             br(),
                             downloadButton('downloadOptMOAO1All', label = "Results",  class="btn-primary")
                           )),
                    column(width = 9,
                           box(
                             title="H1= Heterozygote deficiency: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = "primary",
                             tabsetPanel(type = "tabs",
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Locus", verbatimTextOutput("OptMOAO1outLoc")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Pops", verbatimTextOutput("OptMOAO1outPop"))

                             )
                           ))

)

tabOptMOAO2 = fluidRow(align="left",
                    box(
                      title="Help",
                      tags$iframe(style="height:600px; width:100%",frameborder="0",scrolling="no", src="http://kimura.univ-montp2.fr/%7Erousset/Genepop.pdf#subsubsection.4.1.3"),
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
                             status = "primary",
                             sliderInput("DememoMOAO2",label="Dememorization number", min=100, max = 100000, step=10, value = 10000),
                             sliderInput("NbatchesMOAO2", label="Number of batches", min = 10, max = 10000, step=10, value = 100),
                             sliderInput("NitersMOAO2", label="Number of iterations per batch", min = 400, max = 100000, step=10, value = 5000),
                             actionButton("RunOptMOAO2", label = "Run", icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             status = "primary",
                             style='height:235px;',
                             br(),
                             downloadButton('downloadOptMOAO2All', label = "Results",  class="btn-primary")
                           )),
                    column(width = 9,
                           box(
                             title="H1= Heterozygote excess: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = "primary",
                             tabsetPanel(type = "tabs",
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Locus", verbatimTextOutput("OptMOAO2outLoc")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Pops", verbatimTextOutput("OptMOAO2outPop"))
                             )
                           ))
)

tabOptMOAO3 = fluidRow(align="left",
                    box(
                      title="Help",
                      tags$iframe(style="height:600px; width:100%",frameborder="0",scrolling="no", src="http://kimura.univ-montp2.fr/%7Erousset/Genepop.pdf#subsubsection.4.1.3"),
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
                             status = "primary",
                             sliderInput("DememoMOAO3",label="Dememorization number", min=100, max = 100000, step=10, value = 10000),
                             sliderInput("NbatchesMOAO3", label="Number of batches", min = 10, max = 10000, step=10, value = 100),
                             sliderInput("NitersMOAO3", label="Number of iterations per batch", min = 400, max = 100000, step=10, value = 5000),
                             actionButton("RunOptMOAO3", label = "Run", icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                           ),
                           box(
                             title = "Dowload",
                             width = NULL,
                             solidHeader = TRUE,
                             status = "primary",
                             style='height:235px;',
                             br(),
                             downloadButton('downloadOptMOAO3All', label = "Results",  class="btn-primary")
                           )),
                    column(width = 9,
                           box(
                             title="Probability test: Results",
                             width = NULL,
                             solidHeader = TRUE,
                             status = "primary",
                             tabsetPanel(type = "tabs",
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Locus", verbatimTextOutput("OptMOAO3outLoc")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By Pops", verbatimTextOutput("OptMOAO3outPop")),
                                         tabPanel(style='height:600px; overflow-y: scroll', "By All Locus & All Pops", verbatimTextOutput("OptMOAO3outLocPop"))
                             )
                           ))

)
