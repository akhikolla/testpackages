opt81 <- eventReactive(input$RunOpt81, {
    Gprint(MODE_DEBUG, "Opt81\n")
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        methodAllele = input$methodAllele81
        ncIcoverage = input$coverage81
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(nulls(ficin, outputFile = ficout, nullAlleleMethod = methodAllele, CIcoverage = ncIcoverage), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt81out <- renderText({
    opt <- opt81()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nblig = grep("Number of loci detected", fileText)
            fileText <- paste(fileText[(nblig + 2):length(fileText)], collapse = "\n")
            shinyjs::enable("downloadOpt81All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt81All <- downloadHandler(filename = function() {
    paste("result_opt81_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt81()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
