opt12 <- eventReactive(input$RunOpt12, {
    Gprint(MODE_DEBUG, "Opt12\n")
    dem = input$Dememo12
    nbatchs = input$Nbatches12
    niters = input$Niters12
    ficout = tempfile()
    ficin = GenepopFile()$datapath
    out = TRUE
    if (is.null(ficin)) {
        out = FALSE
    } else {
        Gprint(MODE_DEBUG, ficin)
        setRandomSeed(getSeed(input$randomSeed))
        show("spinner")
        tryCatch(test_HW(ficin, which = "excess", outputFile = ficout, enumeration = TRUE, dememorization = dem, batches = nbatchs, 
            iterations = niters), error = function(e) {
            file.create(ficout)
            write(paste("Exeption : ", e$message), file = ficout)
        }, finally = hide("spinner"))
        file.rename("cmdline.txt", "cmdline.old")
    }
    data.frame(file = ficout, output = out)
})

output$Opt12outLoc <- renderText({
    opt <- opt12()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
            nbli = grep("Results by locus", fileText)
            nblig = grep("Results by population", fileText)
            fileText <- paste(fileText[nbli:(nblig - 1)], collapse = "\n")
            shinyjs::enable("downloadOpt12All")
        } else {
            fileText <- readLines(filePath)
        }
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    fileText
})


output$Opt12outPop <- renderText({
    opt <- opt12()
    if (opt$output) {
        filePath <- toString(opt$file)
        if (file.size(filePath) > 300) {
            fileText <- readLines(filePath)
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

output$downloadOpt12Loc <- downloadHandler(filename = function() {
    paste("result_opt12_Loc_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt12()
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

output$downloadOpt12Pop <- downloadHandler(filename = function() {
    paste("result_opt12_Pop_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt12()
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

output$downloadOpt12All <- downloadHandler(filename = function() {
    paste("result_opt12_", Sys.Date(), ".txt", sep = "")
}, content = function(con) {
    opt <- opt12()
    if (opt$output) {
        filePath <- toString(opt$file)
        fileText <- readLines(filePath)
    } else {
        fileText <- "No genepop file found! please upload a file"
    }
    write(fileText, con)
})
