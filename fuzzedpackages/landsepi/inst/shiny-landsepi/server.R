
simul_params <- createSimulParams(outputDir = paste0(ROOT_PATH,"/www/tmp/"))

## Pathogen parameters
simul_params <- setPathogen(simul_params, loadPathogen("rust"))
## Initial conditions
simul_params <- setInoculum(simul_params, 5e-4)

## Outputs
simul_params <- setOutputs(simul_params, list(epid_outputs="audpc", evol_outputs=""
                                              , thres_breakdown = 50000
                                              , GLAnoDis = 1.48315
                                              , audpc100S = 0.38))


# Server
server <- function(input, output, session) {
  # Functions
  ######################################################################################
  # Test if the cultivars proportion sum is 1
  ProportionValidation <- function () {
    if (input$demo == "RO") {
      sum <-
        ((input$prop0 + input$prop1) + (input$prop0 + input$prop2)) / 2
    }
    else {
      sum <- input$prop0 + input$prop1 + input$prop2
    }
    
    shinyjs::enable(id = "generateLandscape")
    shiny::removeUI(selector = "#propError")
    if (sum != 1 ||
        is.na(sum)) {
      shiny::insertUI(
        selector = "#generateLandscape",
        where = "afterEnd",
        ui = tags$div(
          id = "propError",
          class = "alert alert-warning",
          "The cultivars total proportion should be equal to 1 (100%)"
        )
      )
      shinyjs::disable(id = "generateLandscape")
    }
  }

  #############################################################################
  # Observe EVENT
  #############################################################################

  # About
  # Modal Dialog
  observeEvent(input$About, {
    showModal(modalDialog(
                title = paste0("About : Landsepi V", packageVersion("landsepi")),
                easyClose = TRUE,
                size = "l",
                footer = NULL,
                div(HTML("<h1>Landsepi: Landscape Epidemiology and Evolution</h1>
                          <h3> A stochastic, spatially-explicit, demo-genetic model 
                         simulating the spread and evolution of a plant pathogen in a heterogeneous landscape 
                         to assess resistance deployment strategies. It is based on a spatial geometry for describing 
                         the landscape and allocation of different cultivars, a dispersal kernel for the 
                         dissemination of the pathogen, and a SEIR ('Susceptible-Exposed-Infectious-Removed’) 
                         structure with a discrete time step. It provides a useful tool to assess the performance 
                         of a wide range of deployment options with respect to their epidemiological, 
                         evolutionary and economic outcomes.</h3>
                          <h3> Authors:</h3> J-L Gaussen, J. Papaïx, J-F Rey, L. Rimbaud
                          <h3>Package project:</h3><a href='https://CRAN.R-project.org/package=landsepi' target='_blank'> CRAN package</a><br/><a href='https://gitlab.paca.inra.fr/CSIRO-INRA/landsepi' target='_blank'> Source code</a>
                          <br/> License GPL-3
                          <h3> How to cite the package:</h3> <b>Rimbaud L, Papaïx J, Rey J-F (2019).</b> landsepi: Landscape Epidemiology and Evolution. R package version 0.1.0, &lt;URL: https://cran.r-project.org/package=landsepi&gt;.
                          <h3> Full model description:</h3> <b>Rimbaud L, Papaïx J, Rey J-F, Barrett LG and Thrall PH. 2018.</b> Assessing the durability and efficiency of landscape-based strategies to deploy plant resistance to pathogens. PLoS Computational Biology 14(4): e1006067. <a href='https://doi.org/10.1371/journal.pcbi.1006067' target='_blank'>https://doi.org/10.1371/journal.pcbi.1006067</a>
                          <div>
                          <img src='LogoINRAE_Fr_rvb_web.png' alt='INRAE' style='width:50px; margin-left: 10px;' />
                          <img src='logoBIOSP.jpeg' alt='BioSP' style='width:50px; margin-left: 10px;'/>
                          <img src='PATHO_inra_logo.png' alt='Pathologie végétale' style='width:50px; margin-left: 10px;'/>
                          <img src='CSIRO_Logo.png' alt='CSIRO' style='width:40px; margin-left: 10px;'/>
                          </div>
                "))
    ))
  })

  # Inputs / Outputs
  ######################################################################################
  # Landscape
  shiny::observeEvent(input$landscape, {
    ProportionValidation()
    shinyjs::disable(id = "runSimulation")
  })
  shiny::observeEvent(input$aggregLevel, {
    ProportionValidation()
    shinyjs::disable(id = "runSimulation")
  })

  # Proportion validation
  shiny::observeEvent(input$prop0, {
    ProportionValidation()
    shinyjs::disable(id = "runSimulation")
  })
  shiny::observeEvent(input$prop1, {
    ProportionValidation()
    shinyjs::disable(id = "runSimulation")
  })
  shiny::observeEvent(input$prop2, {
    ProportionValidation()
    shinyjs::disable(id = "runSimulation")
  })
  ######################################################################################
  # Rotation period validation
  shiny::observeEvent(input$rotationPeriod, {
    shinyjs::enable(id = "generateLandscape")
    shinyjs::disable(id = "runSimulation")
    shiny::removeUI(selector = "#rotationPeriodError")
    if (input$demo == "RO") {
      if (input$rotationPeriod < 1 ||
          input$rotationPeriod >= input$nYear ||
          is.na(input$rotationPeriod)) {
        shiny::insertUI(
          selector = "#generateLandscape",
          where = "afterEnd",
          ui = tags$div(
            id = "rotationPeriodError",
            class = "alert alert-warning",
            paste0(
              "The rotation period should be between 1 and ", input$nYear," (the simulation duration)")
          )
        )
        shinyjs::disable(id = "generateLandscape")
      }
    }
  })
  ######################################################################################
  # nYear validation
  shiny::observeEvent(input$nYear, {
    #shinyjs::enable(id = "generateLandscape")
    #shinyjs::enable(id = "runSimulation")
    shinyjs::disable(id = "runSimulation")
    
    shiny::removeUI(selector = "#nYearError")
    if (input$nYear < 1 || input$nYear > 100 || is.na(input$nYear)) {
      shiny::insertUI(
        selector = "#generateLandscape",
        where = "afterEnd",
        ui = tags$div(
          id = "nYearError",
          class = "alert alert-warning",
          paste0(
            "The simulation duration should be between 1 and 100")
        )
      )
      shinyjs::disable(id = "generateLandscape")
      shinyjs::disable(id = "runSimulation")
    }
    else {
      simul_params <<- setTime(simul_params, Nyears = input$nYear, nTSpY = input$nTSpY)
      shinyjs::enable(id = "generateLandscape")}
  })
  ######################################################################################
  # nTSpY validation
  shiny::observeEvent(input$nTSpY, {
    #shinyjs::enable(id = "runSimulation")
    shinyjs::disable(id = "runSimulation")
    
    shiny::removeUI(selector = "#nTSpYError")
    if (input$nTSpY < 1 || input$nTSpY > 365 || is.na(input$nTSpY)) {
      shiny::insertUI(
        selector = "#runSimulation",
        where = "afterEnd",
        ui = tags$div(
          id = "nTSpYError",
          class = "alert alert-warning",
          paste0(
            "The time step should be between 1 and 365")
        )
      )
      shinyjs::disable(id = "runSimulation")
    }
    else {
      simul_params <<- setTime(simul_params, Nyears = input$nYear, nTSpY = input$nTSpY)
      shinyjs::enable(id = "generateLandscape")
    }
  })
  ######################################################################################
  # seed validation
  shiny::observeEvent(input$seed, {
    #shinyjs::enable(id = "generateLandscape")
    #shinyjs::enable(id = "runSimulation")
    shinyjs::disable(id = "runSimulation")
    
    shiny::removeUI(selector = "#seedError")
    if (input$seed < 0 || input$seed > 99999 || is.na(input$seed)) {
      shiny::insertUI(
        selector = "#generateLandscape",
        where = "afterEnd",
        ui = tags$div(
          id = "seedError",
          class = "alert alert-warning",
          paste0(
            "The seed should be between 0 and 99999")
        )
      )
      shinyjs::disable(id = "generateLandscape")
      shinyjs::disable(id = "runSimulation")
    }
    else simul_params <<- setSeed(simul_params, input$seed)
  })
  ######################################################################################
  # Handle the download gpkg button
  output$export <-
    shiny::downloadHandler(filename = "landsepi_landscape.gpkg",
                           content <- function(file) {
                             simul_params <<- saveDeploymentStrategy(simul_params)
                             file.copy(file.path(simul_params@OutputDir, simul_params@OutputGPKG), file)
                           },
                           contentType = "application/x-sqlite3")
  ######################################################################################
  # Handle the "Generate the landscape" button
  shiny::observeEvent(input$generateLandscape, {
    withProgress(message = "Generating Landscape, Please wait...", value = 0, {
      shinyjs::disable(id = "generateLandscape")
      shinyjs::disable(id = "export")
      output$video <- NULL
      #if(!dir.exists(paste0(ROOT_PATH,"/www/tmp/"))) dir.create(paste0(ROOT_PATH,"/www/tmp/"))
      #setwd(paste0(ROOT_PATH,"/www/tmp/"))
      
      # Remove old files
      cleanDir(simul_params@OutputDir)
      
      switch(
        input$demo,
        MO = {
          rotation_period = 0
          rotation_sequence = list(c(simul_params@Croptypes$croptypeID))
          prop = list(c(input$prop0, input$prop1, input$prop2))
          # aggregLevel = strtoi(input$aggregLevel)
        },
        MI = {
          rotation_period = 0
          rotation_sequence = simul_params@Croptypes$croptypeID
          prop = list(c(input$prop0, input$prop1))
          # aggregLevel = strtoi(input$aggregLevel)
        },
        RO = {
          rotation_period = input$rotationPeriod
          prop = list(c(input$prop0, input$prop1),
                      c(input$prop0, input$prop2))
          # aggregLevel = strtoi(input$aggregLevel)
          rotation_sequence <- list(c(simul_params@Croptypes$croptypeID[1], simul_params@Croptypes$croptypeID[2])
                                    , c(simul_params@Croptypes$croptypeID[1], simul_params@Croptypes$croptypeID[3]))
        },
        PY = {
          rotation_sequence = simul_params@Croptypes$croptypeID
          rotation_period = 0
          prop = list(c(input$prop0, input$prop1))
          # aggregLevel = strtoi(input$aggregLevel)
        },
        {
          # Default case
          print("input$generateLandscape : Unknown input$demo")
        }
      )
      simul_params <<- setSeed(simul_params, input$seed)

      incProgress(0.4)
      # Run the landscape generation
      simul_params <<- setLandscape(simul_params, loadLandscape(input$landscape))
      ## Dispersal parameters
      simul_params <<- setDispersalPathogen(simul_params, loadDispersalPathogen(input$landscape))
      ## Dispersal parameters
      disp_host <- loadDispersalHost(simul_params, type = "no")
      simul_params <<- setDispersalHost(simul_params, disp_host)
      
      ## Define the value of aggreg from aggregLevel

      switch(input$aggregLevel
             , "low" = {
               aggreg = 0.07
               algo = "periodic"
               }
             , "medium" = {
               aggreg = 0.25
               algo = "exp"
               }
             , "high" = {
               aggreg = 10
               algo = "periodic"
             }
             , {
               aggreg = 0.25
               algo = "exp"
             }
      )

      simul_params <<- allocateLandscapeCroptypes(simul_params
                                               , rotation_period = rotation_period
                                               , rotation_sequence = rotation_sequence
                                               , rotation_realloc = FALSE
                                               , prop = prop
                                               , aggreg = aggreg
                                               , algo = algo
                                               , graphic = TRUE
      )
      
      setwd(ROOT_PATH)
      
      incProgress(0.5)
      # Print the image of the landscape
      #TODO Loop images for rotation demo
      #output$landscape <- shiny::renderImage({
      #  list(
      #    src = file.path("www/tmp", "landscape_year1.png"),
      #    contentType = 'image/png',
      #    width = "70%",
      #    height = "auto",
      #   alt = "Landscape"
      #  )
      #}, deleteFile = FALSE)
      
      shinyjs::show(id = "landscapeimg")
      output$landscapeimg <- renderPlot({
        imgs <- normalizePath(list.files(simul_params@OutputDir,pattern=".png",full.names = TRUE))
        pngs = lapply(imgs, readPNG)
        asGrobs = lapply(pngs, rasterGrob)
        p <- grid.arrange(grobs=asGrobs, nrow=1 )
        })

      # Using slick : trouble with images size...
      #output$landscape <- renderSlickR({
      #  imgs <- list.files("www",pattern=".png",full.names = TRUE)
      #  slickR(imgs, slickOpts=list(adaptiveHeight=FALSE, respondTo="min"), height = "auto")
      #})
      
      shinyjs::enable(id = "generateLandscape")
      shinyjs::enable(id = "runSimulation")
    })
  })
  
  #############################################################################
  # Stop simulation : Kill future promise
  future_process <- NULL
  observeEvent(input$stopSimulation, {
                 cat(file=stderr(), "STOP button -> stop process ",future_process$job$pid, "\n")
               tools::pskill(future_process$job$pid,signal = tools::SIGTERM)
               tools::pskill(future_process$job$pid,signal = tools::SIGKILL)
  })
 
  ######################################################################################
  # Handle the "Run simulation" button

  shiny::observeEvent(input$runSimulation, {
    withProgress(message = "Running Simulation, please wait...", value = 0 , {
      progressBar <- Progress$new()
      progressBar$set(value = NULL, message = "Running Simulation, please wait...")
      #setwd(paste0(ROOT_PATH,"/www/tmp/"))
      
      shinyjs::disable(id = "inputpanel")
      shinyjs::disable(id = "generateLandscape")
      shinyjs::disable(id = "runSimulation")
      shinyjs::disable(id = "export")
      shinyjs::enable(id = "stopSimulation")
      
      progressBar$set(value = 0.4)

      plan(list(multicore, multiprocess))

      #progressBarTime <- Progress$new()
      #future({
      #  progressBarTime$set(value = 0, message = "please wait...")
      #  for(time in 1:(nYear*1.5)) {
      #    cat(file=stderr(), time/(nYear*1.5), "\n")
      #    progressBarTime$set(value = (time/(nYear*1.5)))
      #    Sys.sleep(1)
      # }
      #  ~progressBarTime$close()
      #
      #})
      #cat(file=stderr(), simul_params)
      #print(simul_params)
      
      future_process <<- future({
        res <- landsepi::runSimul(simul_params,
                                     graphic = FALSE, videoMP4 = TRUE)
      })
      
      then(future_process, onFulfilled = function(value) {
        progressBar$set(value = 0.8, message = "Simulation ended : making video...")
        
        shinyjs::enable(id = "inputpanel")
        shinyjs::enable(id = "generateLandscape")
        #shinyjs::enable(id = "runSimulation")
        shinyjs::enable(id = "export")
        shinyjs::disable(id = "runSimulation")
        shinyjs::disable(id = "stopSimulation")

        output$landscapeimg <- NULL
        hide(id = "landscapeimg")
        
        output$video <-
          shiny::renderUI(
            tags$video(
              id = "video",
              type = "video/mp4",
              src = paste0("tmp/",basename(simul_params@OutputDir),"/video.mp4?rand=",as.integer(Sys.time())),
              controls = "controls",
              width = "100%",
              height = "auto"
            )
          )
        #setwd(dirname(getwd()))
      },
      onRejected = function(err) {
        setwd(ROOT_PATH)
        cat(file=stderr(), "\n ### KILL'EM ALL ### -> Kill simulation \n")
        cleanDir(simul_params@OutputDir)
        shinyjs::enable(id = "inputpanel")
        shinyjs::enable(id = "generateLandscape")
        shinyjs::disable(id = "runSimulation")
        shinyjs::disable(id = "stopSimulation")
        future_process <- NULL
        0
      }) %>%
      finally(~progressBar$close())

     setwd(ROOT_PATH)
    })
  })
  ######################################################################################
  # Handle the demo list
  shiny::observeEvent(input$demo, {
    # Cultivar tab
    switch (
      input$demo,
      MO = {simul_params <<- loadDemoMO(simul_params) },
      MI = {simul_params <<- loadDemoMI(simul_params) },
      RO = { simul_params <<- loadDemoRO(simul_params) },
      PY = { simul_params <<- loadDemoPY(simul_params) },
      {
        # Default case
        print("input$demo : Unknown input$demo")
      }
    )

    output$croptypes <- {
      RenderCroptypes(simul_params@Croptypes)
    }
    output$cultivars <- {
      RenderCultivars(simul_params@CultivarsGenes)
    }
    shinyjs::disable(id = "croptypes")
    shinyjs::disable(id = "cultivars")
    
    # Landscape tab
    shiny::updateSelectInput(session, "landscape", selected = 1)
    simul_params <<- setLandscape(simul_params, loadLandscape(1))
    shiny::updateSelectInput(session, "aggregLevel", selected = "low")
    
    # Enable all the conditionnal inputs by default, we disable it later if needed
    shinyjs::disable(id = "prop2")
    shinyjs::disable(id = "rotationPeriod")
    # Set
    shiny::updateNumericInput(session, "prop2", value = 0)
    shiny::updateNumericInput(session, "rotationPeriod", value = 0)
    
    if (input$demo == "MO") {
      shinyjs::enable(id = "prop2")
      shiny::updateNumericInput(session, "prop0", value = 0.33)
      shiny::updateNumericInput(session, "prop1", value = 0.33)
      shiny::updateNumericInput(session, "prop2", value = 0.34)
      shiny::updateSelectInput(session, "aggregLevel", selected = "high")
    }
    else if (input$demo == "MI" || input$demo == "PY") {
      shiny::updateNumericInput(session, "prop0", value = 0.5)
      shiny::updateNumericInput(session, "prop1", value = 0.5)
      shiny::updateSelectInput(session, "aggregLevel", selected = "high")
    }
    else if (input$demo == "RO") {
      shinyjs::enable(id = "prop2")
      shinyjs::enable(id = "rotationPeriod")
      shiny::updateNumericInput(session, "prop0", value = 0.5)
      shiny::updateNumericInput(session, "prop1", value = 0.5)
      shiny::updateNumericInput(session, "prop2", value = 0.5)
      shiny::updateNumericInput(session, "rotationPeriod", value = 2)
      shiny::updateSelectInput(session, "aggregLevel", selected = "medium")
    }
    else if (input$demo == "PY") {
      shiny::updateSelectInput(session, "aggregLevel", selected = "low")
    }
    
    # More parameters tab
    shiny::updateNumericInput(session, "nYear", value = 30)
    shiny::updateNumericInput(session, "nTSpY", value = 120)
    simul_params <<- setTime(simul_params, Nyears = 30, nTSpY = 120)
    shiny::updateNumericInput(session, "seed", value = 12345)
    simul_params <<- setSeed(simul_params, 12345)
    
    # Disable button
    shinyjs::disable(id = "runSimulation")
    shinyjs::disable(id = "stopSimulation")
    shinyjs::disable(id = "export")
    
    # Remove image
    output$landscapeimg <- NULL
  })
}
