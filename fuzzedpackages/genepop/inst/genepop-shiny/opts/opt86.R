opt86 <- eventReactive(input$RunOpt86, {
    Gprint(MODE_DEBUG, "Opt86\n")
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(sample_haploid(ficin, outputFile = ficout), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt86out <- renderText({
    opt <- opt86()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- paste(readLines(filePath), collapse = "\n")
            shinyjs::enable("downloadOpt86All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt86All <- downloadHandler(filename = function() {
    paste("result_opt86_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt86()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
