opt65 <- eventReactive(input$RunOpt65, {
    Gprint(MODE_DEBUG, "Opt65\n")
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        datatype = input$ploidy65
        statistic = input$isolationStatistic65
        geographicScale = input$geographicScale65
        CIcoverage = input$coverage65
        testPoint = as.numeric(input$testPoint65)
        if (is.na(testPoint)) {
            testPoint = 0
        }
        minimalDistance = input$minDistance65
        maximalDistance = input$maxDistance65
        mantelRankTest = FALSE
        if (input$mantelRankTest65 == "True") {
            mantelRankTest = TRUE
        }
        mantelPermutations = input$mantelPermutations65
        setRandomSeed(getSeed(input$randomSeed))
        setMantelSeed(getMSeed(input$mantelSeed65))
        show("spinner")
        tryCatch(ibd(ficin, outputFile = ficout, settingsFile = "", dataType = datatype, statistic, geographicScale, CIcoverage, 
            testPoint, minimalDistance, maximalDistance, mantelPermutations, mantelRankTest), error = function(e) {
            file.create(ficout)
            write(paste("Exception : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        locus_files = list.files(path = ".", pattern = "LOCUS.")
        for (i in 1:length(locus_files)) {
            file.remove(locus_files[i])
        }
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt65out <- renderText({
    opt <- opt65()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- paste(readLines(filePath), collapse = "\n")
            shinyjs::enable("downloadOpt65All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt65All <- downloadHandler(filename = function() {
    paste("result_opt65_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt65()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
