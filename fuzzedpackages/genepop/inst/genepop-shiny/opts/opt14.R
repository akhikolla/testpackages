# Define a reactive expression for the
opt14 <- eventReactive(input$RunOpt14, {
    Gprint(MODE_DEBUG, "Opt14\n")
    dem = input$Dememo14
    nbatchs = input$Nbatches14
    niters = input$Niters14
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
        tryCatch(test_HW(ficin, which = "global deficit", outputFile = ficout, dememorization = dem, batches = nbatchs, iterations = niters), 
            error = function(e) {
                file.create(ficout)
                write(paste("Exeption : ", e$message), file = ficout)
            }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt14outLoc <- renderText({
    opt <- opt14()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("Results by locus", fileText)
            nblig = grep("Result for all locus and all populations", fileText)
            fileText <- paste(fileText[(nbli + 2):nblig - 2], collapse = "\n")
            shinyjs::enable("downloadOpt14All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found or no space left in the server! please upload a file and try later"
    }
    fileText
})

output$Opt14outPop <- renderText({
    opt <- opt14()
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
        fileText <- "No genepop file found or no space left in the server! please upload a file and try later"
    }
    fileText
})

output$Opt14outLocPop <- renderText({
    opt <- opt14()
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
        fileText <- "No genepop file found or no space left in the server! please upload a file and try later"
    }
    fileText
})


output$downloadOpt14Loc <- downloadHandler(filename = function() {
    paste("result_opt14_Loc_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt14()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
        nbli = grep("Results by locus", fileText)
        nblig = grep("Result for all locus and all populations", fileText)
        fileText <- paste(fileText[(nbli + 2):nblig - 2], collapse = "\n")
    } else {
        fileText <- "No genepop file found or no space left in the server! please upload a file and try later"
    }
    write(fileText, con)
})

output$downloadOpt14Pop <- downloadHandler(filename = function() {
    paste("result_opt14_Pop_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt14()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
        nbli = grep("Results by locus", fileText)
        nblig = grep("Results by population", fileText)
        fileText <- paste(fileText[nblig:(nbli - 2)], collapse = "\n")
    } else {
        fileText <- "No genepop file found or no space left in the server! please upload a file and try later"
    }
    write(fileText, con)
})

output$downloadOpt14All <- downloadHandler(filename = function() {
    paste("result_opt14_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt14()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found or no space left in the server! please upload a file and try later"
    }
    write(fileText, con)
})
