opt34 <- eventReactive(input$RunOpt34, {
    Gprint(MODE_DEBUG, "Opt34\n")
    dem = input$Dememo34
    nbatchs = input$Nbatches34
    niters = input$Niters34
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        show("spinner")
        tryCatch(test_diff(ficin, genic = FALSE, pairs = TRUE, outputFile = ficout, dememorization = dem, batches = nbatchs, 
            iterations = niters), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt34out <- renderText({
    opt <- opt34()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nblig = grep("Iterations per batch", fileText)
            fileText <- paste(fileText[(nblig + 2):length(fileText)], collapse = "\n")
            shinyjs::enable("downloadOpt34All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt34All <- downloadHandler(filename = function() {
    paste("result_opt34_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt34()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
