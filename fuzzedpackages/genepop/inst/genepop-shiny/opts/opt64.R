opt64 <- eventReactive(input$RunOpt64, {
    Gprint(MODE_DEBUG, "Opt64\n")
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        datatype = input$ploidy64
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(Fst(ficin, sizes = TRUE, pairs = TRUE, outputFile = ficout, dataType = datatype), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt64out <- renderText({
    opt <- opt64()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nblig = grep("Number of loci detected", fileText)
            fileText <- paste(fileText[(nblig + 2):length(fileText)], collapse = "\n")
            shinyjs::enable("downloadOpt64All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt64All <- downloadHandler(filename = function() {
    paste("result_opt64_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt64()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
