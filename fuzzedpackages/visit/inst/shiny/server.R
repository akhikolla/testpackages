library(devtools);
load_all();
library(xtable);

shinyServer(function(input, output, session) {
    
    source("visit_ui.R");

    output$mainpage <- renderUI({
        tabpanel.all();
    })
    
    
    output$simulationStart <- renderUI ({
        actionButton(
            inputId = "simu",
            label = "Start simulation"
        )
    })
    
    observeEvent(input$close, {
        stopApp()
    });

    p.data <- reactiveValues();
    
    observeEvent(input$scenarioButton, {
        if (is.null(input$scenarioButton)) return(NULL);
        if (input$scenarioButton == 0) return(NULL);
        
        error <- list();
        
        
        # Probability
        if (input$scenarioInput == 'Probability') {

            sm <- array(0, dim = c(input$ndose, 4));
            for (i in 1:input$ndose) {
                for (j in 1:4) {
                    current <- paste0("predictedOcc", i, j);
                    
                    if (is.null(input[[current]]) || is.na(input[[current]])) {
                        sm[i,j] <- 0;
                    } else {
                        sm[i,j] <- input[[current]];
                    }
                }
            }
            
            if (any(sm < 0)) {
                error <- append(error, "Error: Incorrect values for probability. All values should be nonnegative integers.");
            } else {
                for (i in 1:NROW(sm)) {
                    if (sum(sm[i,]) > 0) {
                        sm[i,] <- sm[i,]/sum(sm[i,]);
                    } else {
                        error <- append(error, "Error: Incorrect values for probability. Each level should contain at least one patient.");
                        break;
                    }
                }
            }
            colnames(sm) <- get.shinyConst()$THETA;
            class(sm) <- get.shinyConst()$CLSTRUEPS;

            
        # Probability by Odds Ratio
        } else if (input$scenarioInput == 'Probability by Odds Ratio') {
            sm <- array(0, dim = c(input$ndose, 3));
            for (i in 1:input$ndose) {
                for (j in 1:3) {
                    current <- paste0("predictedProb", i, j);
                    if (j != 3) {
                        if (is.null(input[[current]]) || is.na(input[[current]])) {
                            sm[i,j] <- NA;
                        } else {
                            sm[i,j] <- input[[current]];
                        }
                    } else {
                        if (input$scenarioRho == 'Multiple') {
                            if (is.null(input[[current]]) || is.na(input[[current]])) {
                                sm[i,j] <- NA;
                            } else {
                                sm[i,j] <- input[[current]];
                            }
                        } else if (input$scenarioRho == 'Single') {
                            if (is.null(input$rho) || is.na(input$rho)) {
                                sm[i,j] <- NA;
                            } else {
                                sm[i,j] <- input$rho;
                            }
                        }
                    }
                }
            }
            
            # Check if p is NA
            if(any(is.na(sm[,1:2]))) {
                error <- append(error, "Error: Missing values for probabilities.");
            }
            
            # Check for 0 <= p <= 1
            if(any(sm[,1:2] < 0, na.rm = TRUE) || any(sm[,1:2] > 1, na.rm = TRUE)) {
                error <- append(error, "Error: Incorrect values for probabilities. Please make sure all probabilities values are between 0 and 1, inclusive.");
            }
            
            # Check if rho is NA
            if(any(is.na(sm[,3]))) {
                error <- append(error, "Error: Missing values for rhos.");
            }
            
            # Check for 0 <= rho
            if(any(sm[,3] < 0, na.rm = TRUE)) {
                error <- append(error, "Error: Incorrect values for rhos. Please make sure all rhos are greater than or equal to 0.");
            }

            if (length(error) == 0) {
                sm <- vtScenario(tox = sm[,1], res = sm[,2], rho = sm[,3]);
            }
        }
        
        
        if (length(error) == 0) {
            
            # Store for rmarkdown
            p.data$sm <- sm;
            
            output$scenario <- renderUI({
                fluidRow(
                    column(12,
                        wellPanel(
                            fluidRow(
                                h4("Scenarios"),
                                style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
                            ),
                            fluidRow(
                                column(6,
                                    renderPlot({
                                        plot(sm, draw.curves = 1:2);
                                    }, bg = 'transparent')
                                ),
                                column(6,
                                    renderPlot({
                                        plot(sm, draw.curves = 3:6);
                                    }, bg = 'transparent')
                                )
                            )
                        )
                    )
                )
            })
            
        } else {
            
            output$scenario <- renderUI({
                return (NULL);
            })
    
            writeError(error)
        }
    })
    
    
    observeEvent(input$simu, {
        if (is.null(input$simu)) return(NULL);
        if (input$simu == 0) return(NULL);
        
        error <- list();
        
        # Check for NA
        decisions <- c(input$dec.cut1, input$dec.cut2, input$dec.cut3, input$etas1, input$etas2)
        if(any(is.na(decisions))) {
            error <- append(error, "Error: Missing values for C1, C2, C3, or DLT boundary values.");
        }
        
        # Check for 0 < d < 1
        if(any(decisions < 0, na.rm = TRUE) || any(decisions > 1, na.rm = TRUE)) {
            error <- append(error, "Error: Incorrect values for C1, C2, C3, or DLT boundary values. Please make sure all values are between 0 and 1, exclusive.");
        }
        
        
        # Unlikely to Happen, but a check for NA probmdl
        if (any(is.na(input$probmdl)) || any(is.null(input$probmdl))) {
            error <- append(error, "Error: No probability model chosen. Please selecte one from Non-Parametric and Non-Parametric+");
        }
        
        
        # Check for PARA and PARA+
        if (input$probmdl == 'PARA' || input$probmdl == 'PARA+') {
            error <- append(error, "Error: Parametric and Parametric+ are unsupported as of this version of shiny website. Please perform the computation through R to bypass this issue.");
        }
        

        # Check for size.level and size.cohort
        if (is.na(input$size.cohort) ||
            is.null(input$size.cohort) ||
            input$size.cohort <= 0 ||
            is.na(input$size.level) ||
            is.null(input$size.level) ||
            input$size.level <= 0) {
            error <- append(error, "Error: Cohort size and level size should be positive integers.");
        } else {
            if (input$size.level <= input$size.cohort) {
                error <- append(error, "Error: Size level should be greater than cohort size.");
            }
        }
        
        
        # Check for n.rep
        if (is.na(input$n.rep) || is.null(input$n.rep) || input$n.rep <= 0) {
            error <- append(error, "Error: Number of replication should be a positive integer.");
        }
        
        # Check for Number of Cores
        n.cores <- input$n.cores;
        if (is.na(input$n.cores) || input$n.cores <= 0 || input$n.cores > parallel::detectCores()) {
            n.cores <- parallel::detectCores() - 1;
        }

    
        # Probability
        if (input$scenarioInput == 'Probability') {
            
            sm <- array(0, dim = c(input$ndose, 4));
            for (i in 1:input$ndose) {
                for (j in 1:4) {
                    current <- paste0("predictedOcc", i, j);
                    
                    if (is.null(input[[current]]) || is.na(input[[current]])) {
                        sm[i,j] <- 0;
                    } else {
                        sm[i,j] <- input[[current]];
                    }
                }
            }
            
            if (any(sm < 0)) {
                error <- append(error, "Error: Incorrect values for probability. All values should be nonnegative integers.");
            } else {
                for (i in 1:NROW(sm)) {
                    if (sum(sm[i,]) > 0) {
                        sm[i,] <- sm[i,]/sum(sm[i,]);
                    } else {
                        error <- append(error, "Error: Incorrect values for probability. Each level should contain at least one patient.");
                        break;
                    }
                }
            }
            colnames(sm) <- get.shinyConst()$THETA;
            class(sm) <- get.shinyConst()$CLSTRUEPS;
            
            
        # Probability by Odds Ratio
        } else if (input$scenarioInput == 'Probability by Odds Ratio') {
            sm <- array(0, dim = c(input$ndose, 3));
            for (i in 1:input$ndose) {
                for (j in 1:3) {
                    current <- paste0("predictedProb", i, j);
                    if (j != 3) {
                        if (is.null(input[[current]]) || is.na(input[[current]])) {
                            sm[i,j] <- NA;
                        } else {
                            sm[i,j] <- input[[current]];
                        }
                    } else {
                        if (input$scenarioRho == 'Multiple') {
                            if (is.null(input[[current]]) || is.na(input[[current]])) {
                                sm[i,j] <- NA;
                            } else {
                                sm[i,j] <- input[[current]];
                            }
                        } else if (input$scenarioRho == 'Single') {
                            if (is.null(input$rho) || is.na(input$rho)) {
                                sm[i,j] <- NA;
                            } else {
                                sm[i,j] <- input$rho;
                            }
                        }
                    }
                }
            }
            
            
            # Check if p is NA
            if(any(is.na(sm[,1:2]))) {
                error <- append(error, "Error: Missing values for probabilities.");
            }
            
            # Check for 0 <= p <= 1
            if(any(sm[,1:2] < 0, na.rm = TRUE) || any(sm[,1:2] > 1, na.rm = TRUE)) {
                error <- append(error, "Error: Incorrect values for probabilities. Please make sure all probabilities values are between 0 and 1, inclusive.");
            }
            
            # Check if rho is NA
            if(any(is.na(sm[,3]))) {
                error <- append(error, "Error: Missing values for rhos.");
            }
            
            # Check for 0 <= rho
            if(any(sm[,3] < 0, na.rm = TRUE)) {
                error <- append(error, "Error: Incorrect values for rhos. Please make sure all rhos are greater than or equal to 0.");
            }
            
            if (length(error) == 0) {
                sm <- vtScenario(tox = sm[,1], res = sm[,2], rho = sm[,3]);
            }
        }
        
        
        if (length(error) == 0) {
            
            progress <- shiny::Progress$new(session, min = 0, max = 1, style = 'old');
            progress$set(message = "Progress", value = 0);
            on.exit(progress$close());

            print(paste0("dec.cut: ", decisions[1:3]));
            print(paste0("etas: ", decisions[4:5]));
            print(paste0("probmdl: ", input$probmdl));
            print(paste0("size.cohort: ", input$size.cohort));
            print(paste0("size.level: ", input$size.level));
            print(paste0("n.rep: ", input$n.rep));
            print(paste0("n.core: ", n.cores));
            print(paste0("seed: ", input$seed1))
            print(paste0("trueps:"));
            print(sm);
            
            showModal(
                modalDialog(
                    id = "prog",
                    title = NULL,
                    size = 's',
                    footer = NULL,
                    fade = FALSE
                )
            )
            
            rst <- vtSimu(
                n.rep = input$n.rep,
                seed = input$seed1,
                trueps = sm,
                size.cohort = input$size.cohort,
                size.level = input$size.level,
                etas = decisions[4:5],
                dec.cut = decisions[1:3],
                prob.mdl = input$probmdl,
                priors = NULL,
                n.cores = n.cores,
                update.progress = progress
            )
   
            
            # Store for rmarkdown
            p.data$rst <- rst;
            
            p.data$ndose <- input$ndose;
            p.data$dec.cut <- c(input$dec.cut1, input$dec.cut2, input$dec.cut3);
            p.data$etas <- c(input$etas1, input$etas2);
            p.data$probmdl <- input$probmdl;
            p.data$size.cohort <- input$size.cohort;
            p.data$size.level <- input$size.level;
            p.data$n.rep <- input$n.rep;
            p.data$n.cores <- n.cores;
            p.data$seed <- input$seed1
            
            output$simulationResult <- renderUI({
                fluidRow(
                    fluidRow(
                        h4("Simulation Result"),
                        style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px; text-align: left;'
                    ),
                    fluidRow(
                        page.simu_output(length(summary(rst)))
                    )
                )
            })
            
            output$simulationStart <- renderUI({
                fluidRow(
                    actionButton(
                        inputId = "simu",
                        label = "Restart simulation"
                    ),
                    align = 'left',
                    style = 'margin-left: 25px; margin-top: 20px !important;'
                )
            })
            
            title <- c("dose: Frequency for each dose level being selected as the optimal dose level",
                "npat: Average number of patients for each cohort and each dose level",
                "samples: Average number of DLT risks and responses for each cohort on each dose level",
                "decision: Frequency each region in the decision map is selected for each cohort on each dose level",
                "prob: Average conditional probabilities corresponding to each region in the decision map for each cohort on each dose level",
                "ptox: Mean and credible interval of DLT risk rates for each cohort on each dose level",
                "pres: Mean and credible interval of immune response rates for each cohort on each dose level");
                
            for (i in 1:(NROW(rst))) {
                local({
                    s <- i;
                    table.name <- paste0("rst.", s);
                    output[[table.name]] <- renderUI ({
                        
                        fluidRow(
                            h3(
                                title[s],
                                style = 'text-align: left;'
                            ),
                            renderTable({
                                xtable(summary(rst)[[s]]);
                            })
                        )
                    })
                })
            }
            removeModal()
            
        } else {
            writeError(error)
        }
    })

    
    observeEvent(input$realButton, {
        if (is.null(input$realButton)) return(NULL);
        if (input$realButton == 0) return(NULL);
        
        error <- list();
        
        
        # Check for NA
        decisions <- c(input$dec.cut1, input$dec.cut2, input$dec.cut3, input$etas1, input$etas2)
        if(any(is.na(decisions))) {
            error <- append(error, "Error: Missing values for C1, C2, C3, or DLT boundary values.");
        }
        
        # Check for 0 < d < 1
        if(any(decisions < 0, na.rm = TRUE) || any(decisions > 1, na.rm = TRUE)) {
            error <- append(error, "Error: Incorrect values for C1, C2, C3, or DLT boundary values. Please make sure all values are between 0 and 1, exclusive.");
        }
        
        
        # Get the matrix m
        om <- array(0, dim = c(input$ndose, 5));
        for (i in 1:input$ndose) {
            om[i,1] <- i;
            for (j in 2:5) {
                current <- paste0("realData", i, j-1);
                if (is.null(input[[current]]) || is.na(input[[current]])) {
                    om[i,j] <- 0;
                } else {
                    om[i,j] <- input[[current]];
                }
            }
        }
        
        # Check for negative numbers in m
        if (any(om < 0)) {
            error <- append(error, "Error: Incorrect values for observed data. Values should be positive integers. ");
        }
        
        
        # Check if input currentLevel is NA
        if (is.na(input$currentLevel) || is.null(input$currentLevel) || input$currentLevel > input$ndose) {
            currentLevel <- input$ndose;
        } else {
            currentLevel <- input$currentLevel;
        }
        
        if (length(error) == 0) {
            
            # Store for rmarkdown
            p.data$om <- om;
            p.data$currentLevel <- currentLevel;
            
            # Generate plotTrack
            output$plotTrack <- renderUI({
                fluidRow(
                    column(8,
                        wellPanel(
                            fluidRow(
                                h4("Track Plot"),
                                style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
                            ),
                            fluidRow(
                                renderPlot({
                                    vtTrack(om, end.width = 0.8, max.level = currentLevel);
                                }, bg = 'transparent')
                            )
                        ),
                        offset = 2
                     )
                 )
            })
            
            # Generate plotInterims
            for (i in 1:(NROW(om))) {
                local({
                    s <- i;
                    plot.name <- paste0("plotInterim", s);
                    
                    output[[plot.name]] <- renderPlot({
                        
                        cur.obs.y <- om[s,-1];
                        if (1 == s) {
                            prev.obs.y <- NULL;
                            prev.res   <- 0;
                        } else {
                            prev.obs.y <- om[s-1,-1];
                            prev.res   <- NULL;
                        }
                        
                        plot(vtInterim(cur.obs.y, prev.obs.y = prev.obs.y, prev.res = prev.res, dec.cut = decisions[1:3], etas = decisions[4:5], seed = input$seed2));
                    }, height = 400, width = 400, bg = 'transparent')
                })
            }
            
            output$plotInterim <- renderUI({
                fluidRow(
                    column(8,
                        wellPanel(
                            fluidRow(
                                h4("Decision Maps"),
                                style = 'margin-left: 20px; border-bottom: 2px solid #E3E3E3; margin-right: 20px;'
                            ),
                            fluidRow(
                                ui.plotInterim(NROW(om))
                            )
                        ),
                        offset = 2
                    )
                )
            })
            
            
        } else {
            
            
            output$plotTrack <- renderUI({
                return (NULL);
            })
            
            for (i in 1:(NROW(om))) {
                local({
                    s <- i;
                    plot.name <- paste0("plotInterim", s);
                    
                    output[[plot.name]] <- renderPlot({
                        return (NULL);
                    })
                })
            }
            
            output$plotInterim <- renderUI({
                return (NULL);
            })
            
            writeError(error)
        }
    }) 

    
    output$export <- downloadHandler(
        filename = function() {
            paste0('export-',
                format(Sys.time(), "%m%d%Y%H%M%S"),
                '.RData'
            )
        },
        content = function(file) {
            export <- get.export();
            save(export, file = file);
        }
    )

    
    get.export <- reactive({
        
        export <- list();
        
        export$ndose <- input$ndose;
        export$dec.cut <- c(input$dec.cut1, input$dec.cut2, input$dec.cut3);
        export$etas <- c(input$etas1, input$etas2);
        export$probmdl <- input$probmdl;
        export$size.cohort <- input$size.cohort;
        export$size.level <- input$size.level;
        export$n.rep <- input$n.rep;
        export$n.cores <- input$n.cores;
        export$currentLevel <- input$currentLevel;
        export$scenarioInput <- input$scenarioInput;
        export$scenarioRho <- input$scenarioRho;
        export$rho <- input$rho;
        export$seed1 <- input$seed1;
        export$seed2 <- input$seed2;
        
        osm <- array(0, dim = c(10, 4));
        om <- array(0, dim = c(10, 4));
        
        for (i in 1:10) {
            for (j in 1:4) {
                current <- paste0("predictedOcc", i, j);
                osm[i,j] <- input[[current]];
                
                current <- paste0("realData", i, j);
                om[i,j] <- input[[current]];
            }
        }
        export$osm <- osm;
        export$om <- om;
        
        psm <- array(0, dim = c(10, 3));
        for (i in 1:10) {
            for (j in 1:3) {
                current <- paste0("predictedProb", i, j);
                psm[i,j] <- input[[current]];
            }
        }
        export$psm <- psm;

        class(export) <- 'visit.export';
        
        return(export);
    })
    
    observeEvent(input$import, {
        if (is.null(input$import)) return (NULL);

        
        if (endsWith(input$import$name, "RData")) {
            load(input$import$datapath)
            if (exists('export') && class(export) == 'visit.export') {
                
                updateSliderInput(
                    session,
                    inputId = "ndose",
                    label = "Number of doses",
                    min = 1,
                    max = 10,
                    value = export$ndose,
                    step = 1
                )
                
                updateNumericInput(
                    session,
                    inputId = "dec.cut1",
                    label = "C1",
                    value = export$dec.cut[1],
                    min = 0,
                    max = 1,
                    step = 0.05
                )
                
                updateNumericInput(
                    session,
                    inputId = "dec.cut2",
                    label = "C2",
                    value = export$dec.cut[2],
                    min = 0,
                    max = 1,
                    step = 0.05
                )
                
                updateNumericInput(
                    session,
                    inputId = "dec.cut3",
                    label = "C3",
                    value = export$dec.cut[3],
                    min = 0,
                    max = 1,
                    step = 0.05
                )
                
                updateNumericInput(
                    session,
                    inputId = "etas1",
                    label = "Lower boundary of DLT risk",
                    value = export$etas[1],
                    min = 0,
                    max = 1,
                    step = 0.05
                )
                
                updateNumericInput(
                    session,
                    inputId = "etas2",
                    label = "Upper boundary of DLT risk",
                    value = export$etas[2],
                    min = 0,
                    max = 1,
                    step = 0.05
                )
                
                updateRadioButtons(
                    session,
                    inputId = "probmdl",
                    label = "",
                    choices = c("Non-Parametric" = "NONPARA",
                        "Non-Parametric+" = "NONPARA+",
                        "Parametric" = "PARA",
                        "Parametric+" = "PARA+"),
                    selected = export$probmdl
                )
                
                updateNumericInput(
                    session,
                    inputId = "size.cohort",
                    label = "Cohort Size",
                    value = export$size.cohort,
                    min = 0,
                    step = 1
                )
                
                updateNumericInput(
                    session,
                    inputId = "size.level",
                    label = "Level Size",
                    value = export$size.level,
                    min = 0,
                    step = 1
                )
                
                updateNumericInput(
                    session,
                    inputId = "n.rep",
                    label = "Number of Replications",
                    value = export$n.rep,
                    min = 1,
                    step = 1
                )
                
                updateNumericInput(
                    session,
                    inputId = "n.cores",
                    label = "Number of Cores",
                    value = export$n.cores,
                    min = 1,
                    step = 1
                )
                
                updateRadioButtons(
                    session,
                    inputId = "scenarioInput",
                    label = "Type",
                    choices = c("Probability by Odds Ratio", "Probability"),
                    selected = export$scenarioInput
                )
                

                updateRadioButtons(
                    session,
                    inputId = "scenarioRho",
                    label = "",
                    choices = c("Single", "Multiple"),
                    selected =  export$scenarioRho
                )


                updateNumericInput(
                    session,
                    inputId = "rho",
                    label = "",
                    value = export$rho,
                    min = 0,
                    step = 1
                )
                
                updateNumericInput(
                    session,
                    inputId = "currentLevel",
                    label = "Current dose level",
                    value = export$currentLevel,
                    min = 1,
                    max = 10,
                    step = 1
                )
                
                updateNumericInput(
                    session,
                    inputId = "seed1",
                    label = "Seed",
                    value = export$seed1,
                    step = 1
                )
                
                updateNumericInput(
                    session,
                    inputId = "seed2",
                    label = "Seed",
                    value = export$seed2,
                    step = 1
                )
               
                
                for (i in 1:10) {
                    for (j in 1:3) {
                        updateNumericInput(
                            session,
                            inputId = paste0("predictedProb", i, j),
                            label = "",
                            value = export$psm[i,j],
                            min = 0,
                            max = 1,
                            step = 0.05
                        )
                    }
                }
                
                for (i in 1:10) {
                    for (j in 1:4) {
                        updateNumericInput(
                            session,
                            inputId = paste0("predictedOcc", i, j),
                            label = "",
                            value = export$osm[i,j],
                            min = 0,
                            step = 1
                        )
                        
                        updateNumericInput(
                            session,
                            inputId = paste0("realData", i, j),
                            label = "",
                            value = export$om[i,j],
                            min = 0,
                            max = 100,
                            step = 1
                        )
                    }
                }
                
                updateTabsetPanel(session, "mainpanel", selected = "About")
                
            } else {
                writeError("Wrong .RData file. Please upload an file that has been exported from VISIT. ")
            }
        } else {
            writeError("Unexpected file. Please import a .RData file.")
        }
    })
    
    output$downloadButton <- downloadHandler(
        filename = function() {
            paste('report_',
                format(Sys.time(), "%m%d%Y%H%M%S"),
                '.',
                switch(input$format,
                    PDF = 'pdf',
                    HTML = 'html',
                    Word = 'docx'
                ),
                sep = ""
            )
        },

        content = function(file) {
            out <- rmarkdown::render('report/report.Rmd',
                switch(input$format,
                    PDF  = rmarkdown::pdf_document(),
                    HTML = rmarkdown::html_document(),
                    Word = rmarkdown::word_document()
                )
            );
            
            bytes <- readBin(out, "raw", file.info(out)$size);
            writeBin(bytes, file);
        }
    )
    
    
    get.data <- reactive({
        
        result <- list();
        
        result$ndose <- p.data$ndose;
        result$dec.cut <- p.data$dec.cut;
        result$etas <- p.data$etas;
        result$probmdl <- p.data$probmdl;
        result$size.cohort <- p.data$size.cohort;
        result$size.level <- p.data$size.level;
        result$n.rep <- p.data$n.rep;
        result$n.cores <- p.data$n.cores;
        result$om <- p.data$om;
        result$currentLevel <- p.data$currentLevel;
        result$sm <- p.data$sm;
        result$rst <- p.data$rst;
        result$seed1 <- p.data$seed1
        
        return(result);
    })
    
})
