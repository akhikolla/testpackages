opt62 <- eventReactive(input$RunOpt62, {
    Gprint(MODE_DEBUG, "Opt62\n")
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        datatype = input$ploidy62
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(Fst(ficin, pairs = TRUE, outputFile = ficout, dataType = datatype), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt62out <- renderText({
    opt <- opt62()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nblig = grep("Number of loci detected", fileText)
            fileText <- paste(fileText[(nblig + 2):length(fileText)], collapse = "\n")
            shinyjs::enable("downloadOpt62All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt62All <- downloadHandler(filename = function() {
    paste("result_opt62_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt62()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
