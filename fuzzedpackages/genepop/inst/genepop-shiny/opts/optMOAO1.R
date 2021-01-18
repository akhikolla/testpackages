optMOAO1 <- eventReactive(input$RunOptMOAO1, {
    Gprint(MODE_DEBUG, "OptMOAO1\n")
    dem = input$DememoMOAO1
    nbatchs = input$NbatchesMOAO1
    niters = input$NitersMOAO1
    # ficout = tempfile()
    ficin = MOAOFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        HWtable_analysis(ficin, which = "deficit", enumeration = TRUE, dememorization = dem, batches = nbatchs, iterations = niters)
        file.rename("cmdline.txt", "cmdline.old")
    }
    cat("1\n")
    filePath <- toString(ficin)
    fileText <- readLines(filePath)
    cat(fileText)
    cat("2\n")
    data.frame(file = ficin, output = out)
})

output$OptMOAO1outLoc <- renderText({
    opt <- optMOAO1()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
        shinyjs::enable("downloadOptMOAO1All")
    } else {
        fileText <- ""
    }
    fileText
})

output$OptMOAO1outPop <- renderText({
    opt <- optMOAO1()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- ""
    }
    fileText
})

output$downloadOptMOAO1Loc <- downloadHandler(filename = function() {
    paste("result_optMOAO1_Loc_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- optMOAO1()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- ""
    }
    write(fileText, con)
})

output$downloadOptMOAO1Pop <- downloadHandler(filename = function() {
    paste("result_optMOAO1_Pop_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- optMOAO1()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- ""
    }
    write(fileText, con)
})

output$downloadOptMOAO1All <- downloadHandler(filename = function() {
    paste("result_optMOAO1_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- optMOAO1()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- ""
    }
    write(fileText, con)
})

