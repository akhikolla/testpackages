

opt66 <- eventReactive(input$RunOpt66, {
    Gprint(MODE_DEBUG, "Opt66\n")
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        datatype = input$ploidy66
        statistic = input$isolationStatistic66
        geographicScale = input$geographicScale66
        CIcoverage = input$coverage66
        testPoint = as.numeric(input$testPoint66)
        if (is.na(testPoint)) {
            testPoint = 0
        }
        minimalDistance = input$minDistance66
        maximalDistance = input$maxDistance66
        mantelRankTest = FALSE
        if (input$mantelRankTest66 == "True") {
            mantelRankTest = TRUE
        }
        mantelPermutations = input$mantelPermutations66
        setRandomSeed(getSeed(input$randomSeed))
        setMantelSeed(getMSeed(input$mantelSeed66))
        show("spinner")
        tryCatch(ibd(ficin, outputFile = ficout, settingsFile = "", dataType = datatype, statistic, geographicScale, CIcoverage, 
            testPoint, minimalDistance, maximalDistance, mantelPermutations, mantelRankTest), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        # locus_files = list.files(path = '.', pattern = 'LOCUS.') for (i in 1:length(locus_files)){ file.remove(locus_files[i]) }
        # ## now removed by the c++ code
        hide("spinner")
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt66out <- renderText({
    opt <- opt66()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- paste(readLines(filePath), collapse = "\n")
            shinyjs::enable("downloadOpt66All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt66All <- downloadHandler(filename = function() {
    paste("result_opt66_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt66()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
