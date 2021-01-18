opt33 <- eventReactive(input$RunOpt33, {
    Gprint(MODE_DEBUG, "Opt33\n")
    dem = input$Dememo33
    nbatchs = input$Nbatches33
    niters = input$Niters33
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(test_diff(ficin, genic = FALSE, pairs = FALSE, outputFile = ficout, dememorization = dem, batches = nbatchs, 
            iterations = niters), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt33out <- renderText({
    opt <- opt33()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nblig = grep("Iterations per batch", fileText)
            fileText <- paste(fileText[(nblig + 2):length(fileText)], collapse = "\n")
            shinyjs::enable("downloadOpt33All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt33All <- downloadHandler(filename = function() {
    paste("result_opt33_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt33()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
