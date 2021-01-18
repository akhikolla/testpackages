opt33 <- eventReactive(input$RunOpt33, {
    cat("Opt33\n")
    dem = input$Dememo33
    nbatchs = input$Nbatches33
    niters = input$Niters33
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    cat(ficin)
    setRandomSeed(getSeed(input$randomSeed))
    RPDGenotypicAllPopulationDifferentiation(ficin, outputFile = ficout, dememorization = dem, batches = nbatchs, iterations = niters)
    file.rename("cmdline.txt", "cmdline.old")
    ficout
})

output$Opt33out <- renderText({
    filePath <- opt33()
    fileText <- paste(readLines(filePath), collapse = "\n")
    fileText
})
