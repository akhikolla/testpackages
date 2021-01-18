opt11 <- eventReactive(input$RunOpt11, {
    Gprint(MODE_DEBUG, "Opt11\n")
    dem <- input$Dememo11
    nbatchs <- input$Nbatches11
    niters <- input$Niters11
    ficout <- tempfile()
    ficin <- GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(test_HW(ficin, which = "deficit", outputFile = ficout, enumeration = TRUE, dememorization = dem, batches = nbatchs,
            iterations = niters), error = function(e) {
            file.create(ficout)
            write(paste("Error: ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt11outLoc <- renderText({
    opt <- opt11()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("Results by locus", fileText)
            nblig = grep("Results by population", fileText)
            fileText <- paste(fileText[nbli:(nblig - 1)], collapse = "\n")
            shinyjs::enable("downloadOpt11All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$Opt11outPop <- renderText({
    opt <- opt11()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            opt11AllValue = fileText
            nblig = grep("Results by population", fileText)
            fileText <- paste(fileText[nblig:length(fileText)], collapse = "\n")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt11Loc <- downloadHandler(filename = function() {
    paste("result_opt11_Loc_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt11()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
        nbli = grep("Results by locus", fileText)
        nblig = grep("Results by population", fileText)
        fileText <- paste(fileText[nbli:(nblig - 1)], collapse = "\n")
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})

output$downloadOpt11Pop <- downloadHandler(filename = function() {
    paste("result_opt11_Pop_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt11()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
        nblig = grep("Results by population", fileText)
        fileText <- paste(fileText[nblig:length(fileText)], collapse = "\n")
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})

output$downloadOpt11All <- downloadHandler(filename = function() {
    paste("result_opt11_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt11()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
