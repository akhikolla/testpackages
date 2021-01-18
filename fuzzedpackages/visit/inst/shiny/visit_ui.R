# Constants in the package
get.shinyConst <- reactive ({
    rst <- list(REGIONS = c("Toxic", "Ineffective", "Safe,Effective", "Effective,Safety concern"),
        THETA     = c("No DLT, No Response", "No DLT, Response",
        "DLT, No Response", "DLT, Response"),
        CLSPRIOR  = "VTPRIOR",
        CLSPOST   = "VTPOST",
        CLSTRUEPS = "VTTRUEPS",
        CLSSIMU   = "VTSIMU",
        CLSDEC    = "VTDEC"
    );
    
    rst$DENLEGEND <- c("Toxicity Rate", "Response Rate", rst$THETA);
    return(rst);
})


##----------------------------------------------------------------------
##                  MAINPAGE UI
##----------------------------------------------------------------------

tabpanel.all <- function() {
    tabsetPanel(type = "pills",
        id = "mainpanel",
        page.about(),
        page.design(),
        page.simu_options(),
        page.simu_result(),
        page.analysis(),
        page.report()
    )
}

##----------------------------------------------------------------------
##                  DISPLAY ERROR
##----------------------------------------------------------------------

writeError <- function(errors) {
    rst <- '';
    for (i in 1:length(errors)) {
        rst <- paste(rst, '<h6>', errors[[i]], '</h6>');
    }
    showModal(
        modalDialog(
            title = "Error",
            footer = NULL,
            size = 'm',
            easyClose = TRUE,
            (HTML(rst))
        )
    )
}


##----------------------------------------------------------------------
##                  ABOUT UI
##----------------------------------------------------------------------

page.about <- function() {
    tabPanel(
        title = "About",
        fluidRow(
            column(8,
                wellPanel(
                    fluidRow(
                        withMathJax(includeHTML('www/text.HTML')),
                        style = 'padding-left: 30px; padding-right: 30px;'
                    ),
                    style = 'padding-left: 30px;'
                ),
                offset = 2
            )
        )
    )
}

##----------------------------------------------------------------------
##                  ABOUT UI
##----------------------------------------------------------------------
page.design <- function() {
    tabPanel(
        title = "Design Options",
        fluidRow(
            column(8,
                page.param(),
                page.prior(),
                offset = 2
            )
        )
    )
}

# Parameters
page.param <- function() {
    wellPanel(
        fluidRow(
            h4("Design Parameters"),
            style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px; margin-bottom: 20px'
        ),
        fluidRow(
            column(6,
                fluidRow(
                    column(10,
                        sliderInput(
                            inputId = "ndose",
                            label = "Number of doses",
                            min = 1,
                            max = 10,
                            value = 3,
                            step = 1
                        ),
                        offset = 1
                    ),
                    style = 'margin-top: 127px; margin-bottom: 127px;'
                ),
                style = 'border-right: 2px solid #E3E3E3;'
            ),
            column(6,
                fluidRow(
                    column(5,
                        numericInput(
                            inputId = "size.cohort",
                            label = "Cohort Size",
                            value = 5,
                            min = 0,
                            step = 1
                        ),
                        style = 'padding: 20px; margin-left: 24px; padding-bottom: 0px; padding-top: 10px;'
                    ),
                    column(5,
                        numericInput(
                            inputId = "size.level",
                            label = "Level Size",
                            value = 10,
                            min = 0,
                            step = 1
                        ),
                        style = 'padding: 20px; margin-left: 17px; padding-bottom: 0px; padding-top: 10px;'
                    )
                ),
                fluidRow(
                    column(5,
                        numericInput(
                            inputId = "etas1",
                            label = "Lower boundary of DLT risk",
                            value = 0.1,
                            min = 0,
                            max = 1,
                            step = 0.05
                        ),
                        style = 'padding: 20px; margin-left: 24px; padding-bottom: 5px;'
                    ),
                    column(5,
                        numericInput(
                            inputId = "etas2",
                            label = "Upper boundary of DLT risk",
                            value = 0.3,
                            min = 0,
                            max = 1,
                            step = 0.05
                        ),
                        style = 'padding: 20px; margin-left: 17px; padding-bottom: 5px;'
                    ),
                    style = 'text-align: center;'
                ),
                fluidRow(
                    column(5,
                        numericInput(
                            inputId = "dec.cut1",
                            label = "C1",
                            value = 0.65,
                            min = 0,
                            max = 1,
                            step = 0.05
                        ),
                        numericInput(
                            inputId = "dec.cut3",
                            label = "C3",
                            value = 0.65,
                            min = 0,
                            max = 1,
                            step = 0.05
                        ),
                        style = 'padding: 20px; margin-left: 24px; padding-top: 5px;'
                    ),
                    column(5,
                        numericInput(
                            inputId = "dec.cut2",
                            label = "C2",
                            value = 0.65,
                            min = 0,
                            max = 1,
                            step = 0.05
                        ),
                        style = 'padding: 20px; margin-left: 17px;padding-top: 5px;'
                    )
                )
            )
        )
    )
}

page.prior <- function() {
    wellPanel(
        fluidRow(
            h4("Probability Model"),
            style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
        ),
        fluidRow(
            column(4,
                radioButtons(
                    inputId = "probmdl",
                    label = "",
                    choices = c("Non-Parametric" = "NONPARA",
                                "Non-Parametric+" = "NONPARA+",
                                "Parametric" = "PARA",
                                "Parametric+" = "PARA+"),
                    selected = "NONPARA"
                ),
                style = 'margin-left: 40px; margin-bottom:40px; margin-top: 20px;'
            ),
            conditionalPanel(
                condition = "(input.probmdl == \"PARA\") || (input.probmdl == \"PARA+\")",
                column(7,
                    fluidRow(
                        column(6,
                            withMathJax(
                                numericInput(
                                    inputId = "vtheta",
                                    label = "$$\\mbox{Variance } \\theta$$",
                                    value = NULL,
                                    min = 1,
                                    step = 1,
                                    width = '90%'
                                )
                            ),
                            withMathJax(
                                numericInput(
                                    inputId = "sdalpha",
                                    label = "$$\\mbox{Standard Deviation } \\alpha$$",
                                    value = NULL,
                                    min = 1,
                                    step = 1,
                                    width = '90%'
                                )
                            ),
                            style = 'margin-top: 12px;'
                        ),
                        column(5,
                            fluidRow(
                                checkboxInput(
                                    inputId = "priory",
                                    label = "Include a prior.y",
                                    value = FALSE
                                ),
                                style = 'margin-top: 33px;'
                            ),
                            offset = 1
                        )
                    )
                )
            ),
            style = 'margin-top: 40px; margin-bottom: 20px; '
        ),
        conditionalPanel(
            condition = "(input.probmdl == \"PARA\") || (input.probmdl == \"PARA+\")",
            fluidRow(
                style = 'border-top: 2px solid #E3E3E3; margin-left: 20px; margin-right: 20px; margin-bottom: 20px;'
            ),
            fluidRow(
                column(2,
                    lapply(1:10, function(i) {
                        fluidRow(
                            conditionalPanel(
                                condition = paste("input.ndose >=", i),
                                column(12,
                                    h3(
                                        paste("Level", i),
                                        style = 'margin-bottom: 23px; margin-top: 24px; padding-top: 0px;'
                                    )
                                )
                            ),
                            style = 'margin-left: 20px;'
                        )
                    }),
                    style = 'margin-top: 54px; margin-left: 0px;'
                ),
                column(10,
                    fluidRow(
                        column(2,
                            h3(
                                paste0(strsplit(get.shinyConst()$THETA[1], ",")[[1]][1], ","),
                                style = 'padding-top: 0px; margin-bottom: 0px;'
                            ),
                            h3(
                                strsplit(get.shinyConst()$THETA[1], ",")[[1]][2],
                                style = 'padding-top: 0px; margin-bottom: 0px;'
                            ),
                            style = 'text-align: center; margin-top: 10px;'
                        ),
                        column(2,
                            h3(
                                paste0(strsplit(get.shinyConst()$THETA[2], ",")[[1]][1], ","),
                                style = 'padding-top: 0px; margin-bottom: 0px;'
                            ),
                            h3(
                                strsplit(get.shinyConst()$THETA[2], ",")[[1]][2],
                                style = 'padding-top: 0px; margin-bottom: 0px;'
                            ),
                            style = 'text-align: center; margin-top: 10px;'
                        ),
                        column(2,
                            h3(
                                paste0(strsplit(get.shinyConst()$THETA[3], ",")[[1]][1], ","),
                                style = 'padding-top: 0px; margin-bottom: 0px;'
                            ),
                            h3(
                                strsplit(get.shinyConst()$THETA[3], ",")[[1]][2],
                                style = 'padding-top: 0px; margin-bottom: 0px;'
                            ),
                            style = 'text-align: center; margin-top: 10px;'
                        ),
                        column(2,
                            h3(
                                paste0(strsplit(get.shinyConst()$THETA[4], ",")[[1]][1], ","),
                                style = 'padding-top: 0px; margin-bottom: 0px;'
                            ),
                            h3(
                                strsplit(get.shinyConst()$THETA[4], ",")[[1]][2],
                                style = 'padding-top: 0px; margin-bottom: 0px;'
                            ),
                            style = 'text-align: center; margin-top: 10px;'
                        ),
                        column(2,
                            "$$\\rho$$",
                            style = 'text-align: center; margin-top: 30px; margin-bottom: 0px;'
                        ),
                        style = 'height: 54px;'
                    ),
                    lapply(1:10, function(i) {
                        fluidRow(
                            conditionalPanel(
                                condition = paste("input.ndose >=", i),
                                lapply(1:4, function(j) {
                                    column(2,
                                        numericInput(
                                            inputId = paste0("prior", i, j),
                                            label = "",
                                            value = NULL,
                                            min = 0,
                                            step = 1)
                                        )
                                }),
                                conditionalPanel(
                                    condition = "(input.probmdl == \"PARA\") || (input.probmdl == \"PARA+\")",
                                    column(2,
                                        numericInput(
                                            inputId = paste0("prior", i, 5),
                                            label = "",
                                            value = NULL,
                                            min = 0,
                                            step = 1
                                        )
                                    )
                                )
                            )
                        )
                    })
                )
            )
        )
    )
}


page.simu_options <- function() {
     tabPanel(
        title = "Simulation Settings",
        fluidRow(
            column(8,
                wellPanel(
                    fluidRow(
                        h4("Settings"),
                        style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
                    ),
                    fluidRow(
                        column(4,
                            numericInput(
                                inputId = "n.rep",
                                label = "Number of Replications",
                                value = 100,
                                min = 1,
                                step = 1
                            )
                        ),
                        column(4,
                            numericInput(
                                inputId = "n.cores",
                                label = "Number of Cores",
                                value = 1,
                                min = 1,
                                step = 1
                            )
                        ),
                        column(4,
                            numericInput(
                                inputId = "seed1",
                                label = "Seed",
                                value = 10000,
                                step = 1
                            )
                        ),
                        style = 'margin-top: 20px; margin-left: 29px; margin-bottom: 8px; margin-right: 29px;'
                    )
                ),
                wellPanel(
                    fluidRow(
                        withMathJax(HTML("<h4>Specify true probabilities: \\(&theta;_{00}, &theta;_{11}, &theta;_{01} , &theta;_{10}\\)</h4>")),
                        style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
                    ),
                    fluidRow(
                        column(4,
                            radioButtons(
                                inputId = "scenarioInput",
                                label = "Type",
                                choices = c("Probability by Odds Ratio", "Probability"),
                                selected = "Probability by Odds Ratio"
                            ),
                            style = 'margin-left: 110px; margin-right: 20px;',
                            offset = 1
                        ),
                        column(4,
                            conditionalPanel(
                                condition = "input.scenarioInput == \"Probability by Odds Ratio\"",
                                fluidRow(
                                    column(2,
                                        withMathJax(
                                            h3(
                                                "$$\\mbox{Number of } \\rho$$",
                                                style = 'margin-top: 20px;'
                                            )
                                        ),
                                        radioButtons(
                                            inputId = "scenarioRho",
                                            label = "",
                                            choices = c("Single", "Multiple"),
                                            selected = "Single"
                                        ),
                                        offset = 5
                                    )
                                )
                            ),
                            style = 'margin-left: 0px; margin-right: 20px;'
                        ),
                        style = 'margin-bottom: 20px;'
                    ),
                    fluidRow(
                        column(2,
                            lapply(1:10, function(i) {
                                fluidRow(
                                    conditionalPanel(
                                        condition = paste("input.ndose >=", i),
                                        column(12,
                                            h3(
                                                paste("Level", i),
                                                style = 'margin-bottom: 23px; margin-top: 24px; padding-top: 0px;'
                                            )
                                        )
                                    ),
                                    style = 'margin-left: 20px;'
                                )
                            }),
                            style = 'margin-top: 54px;'
                        ),
                        
                        column(9,
                            conditionalPanel(
                                condition = "input.scenarioInput == \"Probability by Odds Ratio\"",
                                fluidRow(
                                    column(4,
                                        h3(
                                            "DLT Risk",
                                            style = 'text-align: center; margin-top: 20px; margin-bottom: 0px;'
                                        )
                                    ),
                                    column(4,
                                        h3(
                                            "Immune Response",
                                            style = 'text-align: center; margin-top: 20px; margin-bottom: 0px;'
                                        )
                                    ),
                                    conditionalPanel(
                                        condition = "input.scenarioRho == \"Multiple\"",
                                        withMathJax(
                                            column(4,
                                                "$$\\rho$$",
                                                style = 'text-align: center; margin-top: 30px; margin-bottom: 0px;'
                                            )
                                        )
                                    ),
                                    style = 'height: 54px;'
                                ),
                                lapply(1:10, function(i) {
                                    fluidRow(
                                        conditionalPanel(
                                            condition = paste("input.ndose >=", i),
                                            lapply(1:2, function(j) {
                                                column(4,
                                                    numericInput(
                                                        inputId = paste0("predictedProb", i, j),
                                                        label = "",
                                                        value = 0.2,
                                                        min = 0,
                                                        max = 1,
                                                        step = 0.05
                                                    )
                                                )
                                            }),
                                            column(4,
                                                conditionalPanel(
                                                    condition = "input.scenarioRho == \"Multiple\"",
                                                    numericInput(
                                                        inputId = paste0("predictedProb", i, 3),
                                                        label = "",
                                                        value = NULL,
                                                        min = 0,
                                                        step = 0.05
                                                    )
                                                )
                                            )
                                        )
                                    )
                                }),
                                fluidRow(
                                    column(4,
                                        conditionalPanel(
                                            condition = "input.scenarioRho == \"Single\" && input.scenarioInput == \"Probability by Odds Ratio\"",
                                            
                                            h3(
                                                withMathJax("$$\\rho$$"),
                                                style = 'text-align: center; margin-top: 10px; margin-bottom: 0px;'
                                            ),
                                            
                                            numericInput(
                                                inputId = "rho",
                                                label = "",
                                                value = 1,
                                                min = 0,
                                                step = 1
                                            )
                                        )
                                    )
                                )
                            ),
                            conditionalPanel(
                                condition = "input.scenarioInput == \"Probability\"",
                                fluidRow(
                                    column(3,
                                        h3(
                                            paste0(strsplit(get.shinyConst()$THETA[1], ",")[[1]][1], ","),
                                            style = 'padding-top: 0px; margin-bottom: 0px;'
                                        ),
                                        h3(
                                            strsplit(get.shinyConst()$THETA[1], ",")[[1]][2],
                                            style = 'padding-top: 0px; margin-bottom: 0px;'
                                        ),
                                        style = 'text-align: center; margin-top: 10px;'
                                    ),
                                    column(3,
                                        h3(
                                            paste0(strsplit(get.shinyConst()$THETA[2], ",")[[1]][1], ","),
                                            style = 'padding-top: 0px; margin-bottom: 0px;'
                                        ),
                                        h3(
                                            strsplit(get.shinyConst()$THETA[2], ",")[[1]][2],
                                            style = 'padding-top: 0px; margin-bottom: 0px;'
                                        ),
                                        style = 'text-align: center; margin-top: 10px;'
                                    ),
                                    column(3,
                                        h3(
                                            paste0(strsplit(get.shinyConst()$THETA[3], ",")[[1]][1], ","),
                                            style = 'padding-top: 0px; margin-bottom: 0px;'
                                        ),
                                        h3(
                                            strsplit(get.shinyConst()$THETA[3], ",")[[1]][2],
                                            style = 'padding-top: 0px; margin-bottom: 0px;'
                                        ),
                                        style = 'text-align: center; margin-top: 10px;'
                                    ),
                                    column(3,
                                        h3(
                                            paste0(strsplit(get.shinyConst()$THETA[4], ",")[[1]][1], ","),
                                            style = 'padding-top: 0px; margin-bottom: 0px;'
                                        ),
                                        h3(
                                            strsplit(get.shinyConst()$THETA[4], ",")[[1]][2],
                                            style = 'padding-top: 0px; margin-bottom: 0px;'
                                        ),
                                        style = 'text-align: center; margin-top: 10px;'
                                    ),
                                    style = 'height: 54px;'
                                ),
                                lapply(1:10, function(i) {
                                    fluidRow(
                                        conditionalPanel(
                                            condition = paste("input.ndose >=", i),
                                            lapply(1:4, function(j) {
                                                column(3,
                                                    numericInput(
                                                        inputId = paste0("predictedOcc", i, j),
                                                        label = "",
                                                        value = NULL,
                                                        min = 0,
                                                        step = 1
                                                    )
                                                )
                                            })
                                        )
                                    )
                                })
                            )
                        )
                    ),
                    fluidRow(
                        column(2,
                            actionButton(
                                inputId = "scenarioButton",
                                label = "Simulation Scenario",
                                width = '150px'
                            ),
                            offset = 2
                        ),
                        style = 'margin-top: 30px; margin-bottom: 10px;'
                    )
                ),
                uiOutput("scenario"),
                offset = 2
            )
        )
    )
}

page.simu_result <- function() {
    tabPanel(
        title = "Simulation Result",
        fluidRow(
            column(8,
                wellPanel(
                    uiOutput("simulationResult"),
                    uiOutput("simulationStart"),
                    align = 'center'
                ),
                offset = 2
            )
        )
    )
}



page.simu_output <- function(l) {
    result <- list();
    for (i in 1:l) {
        result[[i]] <- tabPanel(paste("Result", i), uiOutput(paste0("rst.", i)));
    }
    result <- do.call(navlistPanel, c(id = 'r', "", result, well = FALSE, list(widths = c(3, 7))));
    result;
}



##----------------------------------------------------------------------
##                  REAL DATA ANALYSIS
##----------------------------------------------------------------------

page.analysis <- function() {
    tabPanel(
        title = "Real Data Analysis",
        fluidRow(
            column(8,
                wellPanel(
                    fluidRow(
                        h4("Interim Data Analysis"),
                        style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
                    ),
                    fluidRow(
                        column(4,
                            numericInput(
                                inputId = "currentLevel",
                                label = "Current dose level",
                                value = NULL,
                                min = 1,
                                max = 10,
                                step = 1
                            ),
                            numericInput(
                                inputId = "seed2",
                                label = "Seed",
                                value = 10000,
                                step = 1
                            ),
                            style = 'margin-left: 20px;'
                        ),
                        style = 'margin-top: 25px; margin-bottom: 20px;'
                    ),
                    fluidRow(
                        column(2,
                            lapply(1:10, function(i) {
                                fluidRow(
                                    conditionalPanel(
                                        condition = paste("input.ndose >=", i),
                                        column(12,
                                            h3(
                                                paste("Level", i),
                                                style = 'margin-bottom: 23px; margin-top: 24px; padding-top: 0px;'
                                            )
                                        )
                                    ),
                                    style = 'margin-left: 20px;'
                                )
                            }),
                            style = 'margin-top: 54px;'
                        ),
                    
                        column(9,
                            fluidRow(
                                column(3,
                                    h3(
                                        paste0(strsplit(get.shinyConst()$THETA[1], ",")[[1]][1], ","),
                                        style = 'padding-top: 0px; margin-bottom: 0px;'
                                    ),
                                    h3(
                                        strsplit(get.shinyConst()$THETA[1], ",")[[1]][2],
                                        style = 'padding-top: 0px; margin-bottom: 0px;'
                                    ),
                                    style = 'text-align: center; margin-top: 10px;'
                                ),
                                column(3,
                                    h3(
                                        paste0(strsplit(get.shinyConst()$THETA[2], ",")[[1]][1], ","),
                                        style = 'padding-top: 0px; margin-bottom: 0px;'
                                    ),
                                    h3(
                                        strsplit(get.shinyConst()$THETA[2], ",")[[1]][2],
                                        style = 'padding-top: 0px; margin-bottom: 0px;'
                                    ),
                                    style = 'text-align: center; margin-top: 10px;'
                                ),
                                column(3,
                                    h3(
                                        paste0(strsplit(get.shinyConst()$THETA[3], ",")[[1]][1], ","),
                                        style = 'padding-top: 0px; margin-bottom: 0px;'
                                    ),
                                    h3(
                                        strsplit(get.shinyConst()$THETA[3], ",")[[1]][2],
                                        style = 'padding-top: 0px; margin-bottom: 0px;'
                                    ),
                                    style = 'text-align: center; margin-top: 10px;'
                                ),
                                column(3,
                                    h3(
                                        paste0(strsplit(get.shinyConst()$THETA[4], ",")[[1]][1], ","),
                                        style = 'padding-top: 0px; margin-bottom: 0px;'
                                    ),
                                    h3(
                                        strsplit(get.shinyConst()$THETA[4], ",")[[1]][2],
                                        style = 'padding-top: 0px; margin-bottom: 0px;'
                                    ),
                                    style = 'text-align: center; margin-top: 10px;'
                                ),
                                style = 'height: 54px;'
                            ),
                            lapply(1:10, function(i) {
                                fluidRow(
                                    conditionalPanel(
                                        condition = paste("input.ndose >=", i),
                                        lapply(1:4, function(j) {
                                            column(3,
                                                numericInput(
                                                    inputId = paste0("realData", i, j),
                                                    label = "",
                                                    value = 0,
                                                    min = 0,
                                                    max = 100,
                                                    step = 1
                                                )
                                            )
                                        })
                                    )
                                )
                            })
                        )
                    ),
                    
                    fluidRow(
                        column(2,
                            actionButton(
                                inputId = "realButton",
                                label = "Conduct Analysis",
                                width = '135px'
                            ),
                            offset = 2
                        ),
                        style = 'margin-top: 20px;'
                    )
                ),
                offset = 2
            )
        ),
        uiOutput("plotTrack"),
        uiOutput("plotInterim")
    )
}


ui.plotInterim <- function(l){
    result <- list();
    for (i in 1:l) {
        result[[i]] <- tabPanel(paste("Level", i), plotOutput(paste0("plotInterim", i)));
    }
    result <- do.call(navlistPanel, c(id = "t", "", result, well = FALSE, list(widths = c(2, 9))));
    result;
}


##----------------------------------------------------------------------
##                  REPORT
##----------------------------------------------------------------------
page.report <- function() {
    tabPanel(
        title = "Report",
        fluidRow(
            column(8,
                wellPanel(
                    fluidRow(
                        h4("Export and Import Current Settings"),
                        style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
                    ),
                    fluidRow(
                        column(4,
                            downloadButton(
                                outputId = "export",
                                label = "Export",
                                style = 'margin-top: 20px;'
                            ),
                            offset = 2
                        ),
                        column(4,
                            fileInput(
                                inputId = "import",
                                label = "",
                                buttonLabel = "Import"
                            )
                        )
                    )
                ),
                wellPanel(
                    fluidRow(
                        h4("Download Report Analysis"),
                        style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
                    ),
                    fluidRow(
                        radioButtons(
                            inputId = "format",
                            label = "",
                            choices = c('PDF', 'HTML', 'Word')
                        ),
                        downloadButton(
                            outputId = "downloadButton"
                        ),
                        style = 'margin-left: 30px;'
                    )
                ),
                offset = 2
            )
        )
    )
}



