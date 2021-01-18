
opt53 <- eventReactive(input$RunOpt53, {
    Gprint(MODE_DEBUG, "Opt53\n")
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        datatype = input$ploidy53
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(genedivFis(ficin, sizes = TRUE, outputFile = ficout, dataType = datatype), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt53out <- renderText({
    opt <- opt53()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- paste(readLines(filePath), collapse = "\n")
            shinyjs::enable("downloadOpt53All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt53All <- downloadHandler(filename = function() {
    paste("result_opt53_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt53()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
