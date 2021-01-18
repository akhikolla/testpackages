opt32 <- eventReactive(input$RunOpt32, {
    Gprint(MODE_DEBUG, "Opt32\n")
    dem = input$Dememo32
    nbatchs = input$Nbatches32
    niters = input$Niters32
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(test_diff(ficin, genic = TRUE, pairs = TRUE, outputFile = ficout, dememorization = dem, batches = nbatchs, 
            iterations = niters), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt32out <- renderText({
    opt <- opt32()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nblig = grep("Iterations per batch", fileText)
            fileText <- paste(fileText[(nblig + 2):length(fileText)], collapse = "\n")
            shinyjs::enable("downloadOpt32All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt32All <- downloadHandler(filename = function() {
    paste("result_opt32_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt32()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
