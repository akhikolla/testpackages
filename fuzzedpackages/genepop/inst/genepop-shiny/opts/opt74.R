opt74 <- eventReactive(input$RunOpt74, {
    Gprint(MODE_DEBUG, "Opt74\n")
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(conversion(ficin, format = "Linkdos", outputFile = ficout), error = function(e) {
            file.create(ficout)
            write(paste("Exception : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt74out <- renderText({
    opt <- opt74()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- paste(readLines(filePath), collapse = "\n")
            shinyjs::enable("downloadOpt74All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt74All <- downloadHandler(filename = function() {
    paste("result_opt74_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt74()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
