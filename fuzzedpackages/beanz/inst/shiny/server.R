shinyServer(function(input, output, session) {

    userLog                  <- reactiveValues();
    userLog$uid              <- -1;
    userLog$data             <- NULL;
    userLog$SubgroupValid    <- NULL;

    ##load ui functions
    source("beanz_ui.R", local=TRUE);

    ##----------------------------------------------------------------
    ##------------------------MAIN PAGE-------------------------------
    ##----------------------------------------------------------------
    output$mainpage <- renderUI({
        tab.main();
    })

    ##--------------------------------------
    ##---------exit-------------------------
    ##--------------------------------------
    observeEvent(input$close, {stopApp()});


    ##----------------------------------------------------------------
    ##------------------------DATA UPLOAD-----------------------------
    ##----------------------------------------------------------------
    ##--display data uploaded-----
    output$uiData <- DT::renderDataTable({
        if (input$displaydata) {
            get.data();
        }
    }, rownames=NULL, selection="none", options=list(pageLength=50))

    ##get status of fileupload
    output$fileUploaded <- reactive({
       return(is.null(get.data()));
    })

    outputOptions(output,
                  'fileUploaded',
                  suspendWhenHidden=FALSE);

    ##----------------------------------------------------------------
    ##------------------------SPECIFY SUBGROUP------------------------
    ##----------------------------------------------------------------

    ##----------sub group panel-----------
    output$uiSubgroup <- renderUI({
        if (is.null(get.data())) {
            msg.box("Please upload data first.", "warning");
        } else {
            panel.subgroup();
        }
    })

    ##----------select title----------------
    output$txtSel <- renderText({
        paste("Select Variables for ", input$dataformat, sep="");
    });

    ##------RAW:select covariates----------
    output$uiCov <- renderUI({
        dat <- get.data();
        if (is.null(dat) |
            is.null(input$selresp) |
            is.null(input$seltrt))
            return(NULL);

        vname <- colnames(dat);
        vname <- vname[-which(vname == input$selresp | vname == input$seltrt)];

        if (RESP.TYPE[3] == input$resptype) {
            inx.c <- which(vname == input$selcensor);
            if (length(inx.c) > 0)
                vname <- vname[-inx.c];
        }

        checkboxGroupInput("selcov","",
                           choices  = vname,
                           selected = vname);
    })

    ##-----RAW:select response variable----
    output$uiResp <- renderUI({
        dat <- get.data();
        if (is.null(dat))
            return(NULL);

        vname <- colnames(dat);
        selectInput("selresp","",
                    choices = vname,
                    selected = vname[1]);
    })

    ##-----RAW:select treatment variable----
    output$uiTrt <- renderUI({
        dat <- get.data();
        if (is.null(dat))
            return(NULL);

        vname <- colnames(dat);
        selectInput("seltrt","",
                    choices = vname,
                    selected = vname[min(2, length(vname))]);
    })

    ##-----RAW:select censor variable----
    output$uiCensor <- renderUI({
        dat <- get.data();
        if (is.null(dat))
            return(NULL);

        vname <- colnames(get.data());
        selectInput("selcensor","",
                    choices  = vname,
                    selected = vname[min(3, length(vname))]);
    })


    ##------sub:select covariates----------
    output$uisubCov <- renderUI({
        dat <- get.data();
        if (is.null(dat) |
            is.null(input$selsube) |
            is.null(input$selsubvar))
            return(NULL);

        vname <- colnames(dat);
        vname <- vname[-which(vname == input$selsube | vname == input$selsubvar)];

        checkboxGroupInput("selsubcov","",
                           choices  = vname,
                           selected = vname);
    })


    ##-----subgroup:select trt effect variable----
    output$uisubE <- renderUI({
        if (!is.null(get.data())) {
            dat   <- get.data();
            vname <- colnames(dat);
            selectInput("selsube","",
                        choices = vname,
                        selected = vname[1]);
        }
    })

    ##-----subgroup:select treatment variable----
    output$uisubVar <- renderUI({
        dat <- get.data();
        if (is.null(dat) |
            is.null(input$selsube))
            return(NULL);

        vname <- colnames(dat);
        vname <- vname[-which(vname == input$selsube)];
        selectInput("selsubvar","",
                    choices = vname,
                    selected = vname[1]);
    })

    ##------sub:choose nominal covariates----------
    output$uiNominal <- renderUI({
        covs       <- get.sub.cov();
        cov.levels <- get.sub.cov.levels();
        if (is.null(cov.levels))
            return(NULL);

        m.inx <- which(2 < cov.levels);

        if (0 == length(m.inx))
            return(NULL);

        list(h6("Nominal covariates"),
             checkboxGroupInput("selnomcov","",
                                choices  = covs[m.inx],
                                selected = NULL));
    })


    ##------subgroup validation----------
    output$uiValid <- renderText({
        if (0 == input$btnSub)
            return(NULL);

        dat <- get.data();
        if (is.null(dat))
            return(NULL);

        cov.levels <- get.sub.cov.levels();

        if (DATA.FORMAT[1] == input$dataformat) {
            rst <- chk.summary.subgrp(dat, input$selsube, input$selsubvar, input$selsubcov, cov.levels);
        } else {
            rst <- chk.raw.subgrp(dat, input$resptype, input$seltrt, input$selresp,
                                  input$selcensor, input$selcov, cov.levels);
        }

        if (is.null(rst)) {
            rst <- msg.box("Subgroup specification is valid.", "success");
            userLog$SubgroupValid <- 0;
        } else {
            rst <- msg.box(rst, "error");
            userLog$SubgroupValid <- -1;
        }

        HTML(rst);
    })


    ##-----basic subgroup information------
    output$uiBasicdata <- DT::renderDataTable({
        get.subgrp();
    }, options=list(paging=FALSE), rownames=NULL, selection="none")


    ##-----subgroup data----
    output$uiSubgrp <- renderUI({
        if (0 == input$btnSub)
            return(NULL);

        if (-1 == userLog$SubgroupValid)
            return(NULL);

        rst <- get.subgrp();

        if (is.null(rst)) {
            wellPanel(h4("Subgroups"),
                      msg.box("There is error when computing subgroup treatment effect.
                               Please check subgroup specifications.", "error"));
        } else {
            wellPanel(h4("Subgroups"), DT::dataTableOutput("uiBasicdata"));
        }
    })

    ##----------------------------------------------------------------
    ##------------------------MODEL INSTRUCTION-----------------------
    ##----------------------------------------------------------------

    ##-----model instructions-------------
    output$uiMdlInst <- renderText({
        ##mdls <- input$selmdl;
        ##if (is.null(mdls)) return();
        ##rst <- get.inst.header();
        ##for(i in 1:length(ALL.MODELS)) {
        ##  cur.m   <- which(mdls[i] == ALL.MODELS);
        ##  cur.rst <- get.inst.mdl(i);
        ##  rst     <- paste(rst, cur.rst);
        ##}
        fileName <- FMATHMODEL;
        rst      <- readChar(fileName, file.info(fileName)$size);
        HTML(rst);
    })

    ##----------------------------------------------------------------
    ##------------------------CONFIGURATION---------------------------
    ##----------------------------------------------------------------
    output$subgrpsel <- DT::renderDataTable({
        dat <- get.subgrp();
        if (is.null(dat))
            return(NULL);
        dat[, c("Subgroup", get.sub.cov())];
    }, options=list(paging=FALSE, ordering=FALSE, searching=FALSE), rownames=NULL)


    ##----------------------------------------------------------------
    ##------------------------ANALYSIS RESULT-------------------------
    ##----------------------------------------------------------------

    ##----analysis panel------------------
    output$uiAnalysis <- renderUI({
        ##monitor nominal covariates
        input$selnomcov;

        ##monitor selected models
        get.ana.models();

        if (is.null(get.data())) {
            msg.box("Please upload data first.", "warning");
        } else if (is.null(get.subgrp())) {
            msg.box("Please specify subgroups first.", "warning");
        } else {
            list(msg.warnings(),
                 msg.box("<p>Click button below to conduct Bayesian analysis.</p>"),
                 actionButton("btnAna", "Conduct Bayesian Analysis", styleclass = "success"),
                 uiOutput("uiRhatWarn"),
                 uiOutput("mdlrst"));
        };
    })

    ##check convergence
    output$uiRhatWarn <- renderUI({
        arst <- ana.rst();
        if (is.null(arst))
            return(NULL);
        flag <- FALSE;
        for (i in 1:length(arst)) {
            cur.rhat <- arst[[i]]$rhat;
            inx      <- which(cur.rhat > get.rhat.warn());
            if (0 < length(inx)) {
                flag <- TRUE;
                break;
            }
        }
        if (!flag) {
            return(NULL);
        } else {
            rst <- msg.box("There are models that may have convergence issues.
                  Please check Rhat of each model for
                  details. Longer iterations are recommended to solve this issue. If the issue remains,
                  it is recommended to carefuly check the convergence by techniques such as
                  Brook-Gelman-Rubin plots.", "error");
        }
        rst
    })


    ##---output results panel---
    output$mdlrst <- renderUI({
        if (is.null(input$btnAna))
            return(NULL);

        if (0 == input$btnAna)
            return(NULL);

        if (is.null(ana.rst()))
            return(NULL);

        if (input$displayby == DISPLAY.BY[1]) {
            do.call(navlistPanel, gen.rst.tabs());
        } else {
            do.call(navlistPanel, gen.rst.tabs.mdl());
        };
    })


    ##----present stan results---------------------
    observe({
        if (is.null(input$btnAna))
            return(NULL);

        sub.sel       <- get.subgrp.selected();
        disp.cut      <- input$displaycut;
        disp.cutcomp  <- input$displaycutcomp;
        disp.digit    <- input$displaydigit;
        disp.ref      <- input$displayref;

        isolate({
            arst <- ana.rst();

            if (is.null(arst))
                return(NULL);

            ##predictive distribution
            apred <- pred.rst();
            ##current subgroup data
            dat.sub <- get.subgrp();

            ##models dic
            dic <- get.all.dic();

            for (i in 1:length(arst)) {
                local({
                    myi  <- i;
                    ##add reference
                    if (1 == myi | !disp.ref) {
                        ref <- NULL;
                    } else {
                        ref <- arst[[1]];
                    }

                    ##----effect
                    plotname <- paste("plotrst", myi, sep="");
                    output[[plotname]] <- renderPlot({
                        bzPlot(arst[[myi]], sel.grps=sub.sel, ref.stan.rst=ref);
                    }, bg="transparent", width=900, height=450)

                    plotname <- paste("plotforest", myi, sep="");
                    output[[plotname]] <- renderPlot({
                        bzForest(arst[[myi]], sel.grps=sub.sel, cut=disp.cut, ref.stan.rst=ref);
                    }, bg="transparent", width=900, height=450)

                    plotname <- paste("plotpred", myi, sep="");
                    output[[plotname]] <- renderPlot({
                        plot.pred(apred[[myi]], dat.sub, SUB.HEAD);
                    }, bg="transparent", width=900, height=450)

                    uname <- paste("tblrst", myi, sep="");
                    output[[uname]] <- DT::renderDataTable(
                        {
                            bzSummary(arst[[myi]], sel.grps=sub.sel, cut=disp.cut,
                                           digits=disp.digit, ref.stan.rst=ref);
                        },
                        rownames=NULL,
                        selection="none",
                        options=list(paging=FALSE)
                    );

                    ##----comparison
                    plotname <- paste("plotcomp", myi, sep="");
                    output[[plotname]] <- renderPlot({
                        bzPlotComp(arst[[myi]], sel.grps=sub.sel);
                    }, bg="transparent")

                    plotname <- paste("plotforestcomp", myi, sep="");
                    output[[plotname]] <- renderPlot({
                        bzForestComp(arst[[myi]], sel.grps=sub.sel, cut=disp.cutcomp);
                    }, bg="transparent", width=500, height=400);

                    uname <- paste("tblcomp", myi, sep="");
                    output[[uname]] <- DT::renderDataTable(
                        {
                            bzSummaryComp(arst[[myi]], sel.grps=sub.sel,
                                              cut=disp.cutcomp, digits=disp.digit);
                        },
                        rownames=NULL,
                        selection="none",
                        options=list(paging=FALSE)
                    );

                    ##------stan
                    uname <- paste("txtdic", myi, sep="");
                    output[[uname]] <- renderText({
                        rst <- sprintf("<h6> <p>The  leave-one-out cross-validation information
                                    criterion (LOOIC) of this model is <strong>%5.3f</strong>. </p>
                                    <p> The deviance information criterion(DIC) of
                                    this model is <strong>%5.3f</strong>.</p>
                                   </h6>", dic[myi, 2], dic[myi, 1]);
                        HTML(rst);
                    })

                    uname <- paste("txtrst", myi, sep="");
                    output[[uname]] <- renderPrint({
                        print(arst[[myi]]$stan.rst);
                    })

                    uname <- paste("plottrace", myi, sep="");
                    output[[uname]] <- renderPlot({
                        cur.stan <- arst[[myi]]$stan.rst;
                        rstan::traceplot(cur.stan, pars = get.pars(cur.stan@model_pars));
                        ##r.plot.trace(arst[[myi]]);
                    },
                    bg="transparent",
                    height=get.traceplot.height(arst[[myi]]$smps))

                    uname <- paste("plotrhat", myi, sep="");
                    output[[uname]] <- renderPlot({
                        cur.stan <- arst[[myi]]$stan.rst;
                        rstan::stan_rhat(cur.stan);
                    },
                    bg="transparent")

                    uname <- paste("txtrhat", myi, sep="");
                    output[[uname]] <- renderText({
                        cur.rhat <- arst[[myi]]$rhat;
                        inx      <- which(cur.rhat > get.rhat.warn());
                        if (0 == length(inx)) {
                            ss <- "None";
                        } else {
                            ss <- paste(names(cur.rhat[inx], collapse=","));
                        }
                        rst <- sprintf("<h6><p>Rhat is the Gelman and Rubin potential scale
                                        reduction statistic that
                                        measures whether the sampling chains have converged. When If the chains
                                        have not converged to a common distribution, the Rhat statistic will be
                                        greater than one. </p> <p> The following parameters have Rhat bigger
                                        than 1.1: <strong>%s</strong>.  </p>",
                                       ss);
                        HTML(rst);
                    })

                    uname <- paste("plotrhat", myi, sep="");
                    output[[uname]] <- renderPlot({
                        cur.stan <- arst[[myi]]$stan.rst;
                        rstan::stan_rhat(cur.stan);
                    },
                    bg="transparent")

                    ## gelman-rubin plot
                    uname <- paste("plotgr", myi, sep="");
                    output[[uname]] <- renderPlot({
                        cur.stan  <- arst[[myi]]$stan.rst;
                        cur.array <- as.array(cur.stan);
                        cur.pars  <- get.pars(attr(cur.array, "dimnames")$parameters,
                                              "^mu+");
                        cur.coda  <- coda::mcmc.list(lapply(1:ncol(cur.stan),
                                                           function(x)
                                              coda::mcmc(cur.array[,x,cur.pars])));
                        coda::gelman.plot(cur.coda);
                    },
                    bg="transparent",
                    height=get.traceplot.height(arst[[myi]]$smps))

                })
            }
        })
    })

    ##----------------------------------------------------------------
    ##------------------------ToolBOX---------------------------------
    ##----------------------------------------------------------------

    ##----------toolbox panel-----------
    output$uiTool <- renderUI({
        if (is.null(get.data())) {
            rst <- msg.box("Please upload data first.", "warning");
        } else if (is.null(get.subgrp())) {
            rst <- msg.box("Please specify subgroups first.", "warning");
        } else {
            lele <- list(includeHTML("www/beanz_anoint.html"));
            if (DATA.FORMAT[1] == input$dataformat) {
                lele <- c(lele,
                          list(msg.box("Toolbox function ANOINT only works for
                                        subject level raw data.", "warning")
                               ));
            } else {
                if (!requireNamespace("anoint", quietly = TRUE)) {
                    lele <- c(lele,
                              list(wellPanel(
                                  msg.box("R package anoint is needed for this function
                                           to work. Please install.",
                                          "warning")
                              )));
                } else {
                    lele <- c(lele,
                              list(withTags(
                                  div(actionButton("btnTool", "Conduct Analysis", styleclass = "success"),
                                      style="margin-bottom:15px")),
                                  uiOutput("uiAnoint"))
                              );
                }
            };

            rst <- list(msg.warnings(),
                        do.call(wellPanel, lele));
            rst <- c(rst,
                     list(wellPanel(includeHTML("www/beanz_gailsimon.html"),
                                    fluidRow(
                                        column(3, h6("Clinically meaningful threshold"),
                                               numericInput(inputId="inGscut", label="",
                                                            value=0, min=0, step=0.1))),
                                    htmlOutput("txtGS")
                                    )));
        }

        rst
    })


    output$txtGS <- renderText({
        dat <- get.subgrp();
        if (is.null(dat)) {
            rst  <- paste("<h6>Subgroups are not correctly specified. Please check</h6>");
        } else {
            pval <- bzGailSimon(dat[, SUB.HEAD[1]], sqrt(dat[, SUB.HEAD[2]]), input$inGscut); ;
            rst  <- sprintf("<h6> The P-value for testing the null hypothesis of no qualitative
                                  interaction is %5.4f.</h6>", pval);
        }
        HTML(rst);
    })

    output$uiAnoint <- renderUI({
        if (0 == input$btnTool)
            return(NULL);

        tabsetPanel(type = "pills",
                    tabPanel("One-by-one interaction",
                             verbatimTextOutput("txtOio")),
                    tabPanel("Unstructured interaction",
                             verbatimTextOutput("txtUim")),
                    tabPanel("Proportional interaction",
                             verbatimTextOutput("txtPim"))
                    )
    })

    output$txtOio <- renderPrint({
        ano <- get.anoint();
        if (is.null(ano))
            return(NULL);
        print(ano$obo);
    })

    output$txtUim <- renderPrint({
        ano <- get.anoint();
        if (is.null(ano))
            return(NULL);
        print(ano$uim);
    })

    output$txtPim <- renderPrint({
        ano <- get.anoint();
        if (is.null(ano))
            return(NULL);
        print(ano$pim);
    })


    ##----------------------------------------------------------------
    ##------------------------DOWNLOAD--------------------------------
    ##----------------------------------------------------------------

    ##report tab
    output$uiReport <- renderUI({
        if (is.null(get.data())) {
            msg.box("Please upload data first.", "warning");
        } else if (is.null(get.subgrp())) {
            msg.box("Please specify subgroups first.", "warning");
        } else if (is.null(ana.rst())) {
            msg.box("Please conduct analysis first.", "warning");
        } else {
            panel.download();
        }
    })


    ##summary of the selected model
    output$rstsummary <- renderText({
        dic <- get.all.dic();

        if (is.null(dic))
            return(NULL);

        min.mdl <- which.min(dic[,2]);
        mdls    <- get.ana.models();
        rst <- paste("<h6> Based on leave-one-out cross-validation information
                      criterion (LOOIC), results from the <strong>", mdls[min.mdl], " </strong>
                      model are reported here. <strong>However, this does not imply that the selected
                      model is significantly better than the other models </strong>. Please consider the
                      totality of the information for drawing the conclusions. </h6>");
        HTML(rst);
    })

    ##selected model results
    output$rsttbl <- DT::renderDataTable({
        arst       <- ana.rst();
        dat        <- get.subgrp();
        dic        <- get.all.dic();
        disp.cut   <- input$displaycut;
        disp.digit <- input$displaydigit;
        sub.cov    <- get.sub.cov();
        rst        <- bzRptTbl(arst, dat, sub.cov, disp.cut, disp.digit);
    }, options=list(paging=FALSE), rownames=NULL, selection="none")


    ##----download -----
    output$btnDload <- downloadHandler(
        filename=function() {
            paste('report_',
                  format(Sys.time(), "%m%d%Y%H%M%S"),
                  '.',
                  switch(input$format,
                         PDF = 'pdf',
                         HTML = 'html',
                         Word = 'docx'
                         ),
                  sep="")
        },

        content=function(file) {
            out <- rmarkdown::render('report/report.Rmd',
                           switch(input$format,
                                  PDF  = rmarkdown::pdf_document(),
                                  HTML = rmarkdown::html_document(),
                                  Word = rmarkdown::word_document()
                          ));

            bytes <- readBin(out, "raw", file.info(out)$size);
            writeBin(bytes, file);
        })

})
