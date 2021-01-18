opt15 <- eventReactive(input$RunOpt15, {
    Gprint(MODE_DEBUG, "Opt15\n")
    dem = input$Dememo15
    nbatchs = input$Nbatches15
    niters = input$Niters15
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    du <- disk.usage()/(1024 * 1024)
    # here we nedd at leat 100 Go of free space to allow temporary files to be
    if (is.null(ficin) || du[2] < 100) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(test_HW(ficin, which = "global excess", outputFile = ficout, dememorization = dem, batches = nbatchs, iterations = niters), 
            error = function(e) {
                file.create(ficout)
                write(paste("Exeption : ", e$message), file = ficout)
            }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt15outLoc <- renderText({
    opt <- opt15()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("Results by locus", fileText)
            nblig = grep("Result for all locus and all populations", fileText)
            fileText <- paste(fileText[(nbli + 2):nblig - 2], collapse = "\n")
            shinyjs::enable("downloadOpt15All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$Opt15outPop <- renderText({
    opt <- opt15()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("Results by locus", fileText)
            nblig = grep("Results by population", fileText)
            fileText <- paste(fileText[nblig:(nbli - 2)], collapse = "\n")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$Opt15outLocPop <- renderText({
    opt <- opt15()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("Result for all locus and all populations", fileText)
            fileText <- paste(fileText[nbli:length(fileText)], collapse = "\n")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})

output$downloadOpt15Loc <- downloadHandler(filename = function() {
    paste("result_opt15_Loc_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt15()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
        nbli = grep("Results by locus", fileText)
        nblig = grep("Result for all locus and all populations", fileText)
        fileText <- paste(fileText[(nbli + 2):nblig - 2], collapse = "\n")
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})

output$downloadOpt15Pop <- downloadHandler(filename = function() {
    paste("result_opt15_Pop_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt15()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
        nbli = grep("Results by locus", fileText)
        nblig = grep("Results by population", fileText)
        fileText <- paste(fileText[nblig:(nbli - 2)], collapse = "\n")
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})

output$downloadOpt15All <- downloadHandler(filename = function() {
    paste("result_opt15_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt15()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
