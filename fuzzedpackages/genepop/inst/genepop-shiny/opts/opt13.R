opt13 <- eventReactive(input$RunOpt13, {
    Gprint(MODE_DEBUG, "Opt13\n")
    dem = input$Dememo13
    nbatchs = input$Nbatches13
    niters = input$Niters13
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(test_HW(ficin, which = "Proba", outputFile = ficout, enumeration = TRUE, dememorization = dem, batches = nbatchs, 
            iterations = niters), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt13outLoc <- renderText({
    opt <- opt13()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("Results by locus", fileText)
            nblig = grep("Results by population", fileText)
            fileText <- paste(fileText[nbli:(nblig - 1)], collapse = "\n")
            shinyjs::enable("downloadOpt13All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$Opt13outPop <- renderText({
    opt <- opt13()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("All locus, all populations", fileText)
            nblig = grep("Results by population", fileText)
            fileText <- paste(fileText[nblig:(nbli - 1)], collapse = "\n")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$Opt13outLocPop <- renderText({
    opt <- opt13()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("All locus, all populations", fileText)
            fileText <- paste(fileText[nbli:length(fileText)], collapse = "\n")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt13Loc <- downloadHandler(filename = function() {
    paste("result_opt13_Loc_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt13()
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

output$downloadOpt13Pop <- downloadHandler(filename = function() {
    paste("result_opt13_Pop_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt13()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
        nbli = grep("All locus, all populations", fileText)
        nblig = grep("Results by population", fileText)
        fileText <- paste(fileText[nblig:(nbli - 1)], collapse = "\n")
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})

output$downloadOpt13All <- downloadHandler(filename = function() {
    paste("result_opt13_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt13()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
