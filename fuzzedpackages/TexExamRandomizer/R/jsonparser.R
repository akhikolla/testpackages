# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017
#### Parsing document for JSON TexExamRandomizer directions
####
####

#' @title ParsePreambleForOptions
#'
#' @description This function parses a preamble of a document trying to read options handed to the package TexExamRandomizer to be used in compiling.
#'
#' @param preamble character vector identifying the preamble from which to pass the JSON readon through
#'
#' @return Returns a list, that concatenates all the lists of options described on the file.
#'
#' @details It find all \code{\%!TexExamRandomizer = {}} lines. It then uses the function  \code{\link[jsonlite]{fromJSON}} to parse them, and it concatenates all those options.
#'
#' If more than one option with the same name is given, it tries to concatenate those.
#' However, it doesn't do that recursively, only if the names of the outer layer are the same... therefore, in nested structure you might end up with a list that have twice the same name. Keep in mind that in those cases, the default behaviour of R is to select the first one.
#'
#' @family jsoncompiler
#'
#' @export
ParsePreambleForOptions <- function(preamble) {
    connectList <- function(A, B) {
        #Connects two lists if they have conflicting names
        conflicting_names <- names(B) %in% names(A)
        for (name in names(B)[conflicting_names]) {
            A[[name]] <- c(A[[name]], B[[name]])
        }

        return(c(A, B[!conflicting_names]))
    }


    options_pattern <- "(?i)%\\s*!\\s*TexExamRandomizer\\s*=\\s*(\\{.*\\})"


    json_info <-
        stats::na.omit(
            stringr::str_match(
                preamble,
                options_pattern
            )[,2]
        )
    # First column is the non-matches, second column is the first group, which is greedy... so it should grab everything, assuming someone didn't write something to disturb it


    optionList <- list()

    for (line in json_info) {

         tryCatch(
            {
                optionList <-
                    connectList(
                        optionList,
                        jsonlite::fromJSON(line)
                    )
            },
            error = function(e) {
                cat("JSON couldn't be read from:\t", line, "\nProceeding ignoring this\n")
            }

        )
    }


    return(optionList)

}



#' @title Json Exam Document Parser
#'
#' @description This function takes a series of options as obtained from \code{\link[optparse]{parse_args}} through the parameter \code{opt}. The "examples" section provides all the options that it can parse.
#'
#' From within those options, a \code{--file} option is mandatory.
#'
#' The file option provides a 'LaTeX' file name in which to search for lines on the preamble \code{\%!TexExamRandomizer} within the first 200 lines.
#'
#' With those options that it finds through tags, it passes the function \code{\link{CreateRandomExams}}.
#'
#' Note that the tags must respect the JSON format, that is. It \emph{needs} to be written within double quotes.
#' @param opt Options as parsed from \code{\link[optparse]{parse_args}}. The function expects a series of options, the example code exemplifies those options that the function understands.
#'
#'
#' @details
#'
#' All the options can be found on
#'
#' \code{vignette("ExamOptions", package = "TexExamRandomizer")}
#'
#' The options that are called "command line" options in the vignette are those that are given to the function through \code{opt}, the rest of the options are read directly froo the document specified with \code{--file <filename>}
#'
#'
#' @family jsoncompiler
#'
#' @export
#' @examples
#' \dontrun{
#' #!/bin/Rscript
#' #This example showcases the type of script this jsonparser might be used on.
#' # You can still use it without a script,
#' # just by adding a list that has the same names as the list provided in opt
#' library(optparse)
#' option_list <- list(
#'     make_option(
#'         c("--file"),
#'         action = "store",
#'         default = NULL,
#'         type = 'character',
#'         help = "Filename of the Tex File"
#'     ),
#'     make_option(
#'         c("--table"),
#'         action = "store",
#'         default = NULL,
#'         type = 'character',
#'         help = "Filename of the table to break down. It overwrites the values written on the file"
#'     ),
#'     make_option(
#'         c("-n", "--noutput"),
#'         action = "store",
#'         default = NULL,
#'         type = "integer",
#'         help = "Number of output Versions"
#'     ),
#'     make_option(
#'         c("-q", "--nquestions"),
#'         action = "store",
#'         default = NULL,
#'         type = "character",
#'         help = "Number of output questions"
#'     ),
#'     make_option(
#'         c("-s", "--seed"),
#'         action = "store",
#'         default = NULL,
#'         type = "integer",
#'         help = "Seed for any randomization done"
#'     ),
#'     make_option(
#'         c("-c", "--compile"),
#'         action = "store_true",
#'         default = FALSE,
#'         type = "logical",
#'         help = "Should the output folder be compiled or not"
#'     ),
#'     make_option(
#'         c("--xelatex"),
#'         action = "store_true",
#'         default = FALSE,
#'         type = "logical",
#'         help = "Should we use xelatex to compile or not"
#'     ),
#'     make_option(
#'         c("-d", "--debug"),
#'         action = "store_true",
#'         default = FALSE,
#'         type = "logical",
#'         help = "If debugging, it doesn't remove auxiliary files"
#'     )
#' )
#'
#'
#' #### PARSING OPTIONS ####
#' ####
#' opt <-
#'     parse_args(
#'         OptionParser(option_list = option_list),
#'         positional_arguments = TRUE
#'     )
#'
#' # Let's assume the file was the example file
#' testfile <-
#'     system.file(
#'         "extdata",
#'         "ExampleTexDocuments",
#'         "exam_testing_nquestions.tex", #Test exam that doesn't require a table
#'         package = "TexExamRandomizer")
#'
#' # To prevent modifying the file system in examples
#' temporalfile <- paste(tempfile(), ".tex", sep = "")
#'
#' file.copy(testfile, temporalfile)
#'
#' opt$options$file <- temporalfile
#'
#'
#'
#' jsonexamparser(opt)
#' }
#'

jsonexamparser <- function(opt) {

    # Checking main tex file exists

    if (is.na(opt$options$file)) {
        stop("Main Tex file not passed")
    } else {
        texFile <- opt$options$file
        cat(
            texFile,
            '...exist =',
            assertthat::assert_that(file.exists(texFile)),
            "\n",
            file = stderr()
        )

    }

    MainDirectory <- dirname(texFile)

    compileInfo <- "compileinfo.txt" # This will hold the information from the compilation steps

    fullAnswerSheet <- "fullanswersheet.csv" #Filename with the full answer sheet

    shortAnswerSheet <- "shortanswersheet.csv" #Filename with the short version of the answer sheet

    this.basename <- basename(texFile)
    assertthat::assert_that(assertthat::has_extension(texFile, "tex"))
    outputBaseName <- sub(pattern = "\\.tex$", replacement = "", x = this.basename)
    outputDirectory <- file.path(MainDirectory, outputBaseName)


    x <- readLines(texFile)


    # We only search for the %!TexExamRandomizer within the first 200 lines of the document
    preamble_options <- ParsePreambleForOptions(x[1:200])


    if (is.null(opt$options$table)) {
        tableFile <- preamble_options$table
    } else {
        tableFile <- opt$options$table
    }


    if (is.null(opt$options$noutput)) {
        nOutputVersions <- preamble_options$noutput
    } else {
        nOutputVersions <- opt$options$noutput
    }

    if (is.null(opt$options$noutput)) {
        nOutputQuestions <- preamble_options$nquestions
    } else {
        nOutputQuestions <- opt$options$nquestions
    }



    if (!is.null(opt$options$seed)) {
        seed <- opt$options$seed
    } else if (!is.null(preamble_options$seed)) {
        seed <- preamble_options$seed
    } else if (file.exists(file.path(outputDirectory, compileInfo))) {
        seed <- NULL

        tryCatch(
            {
                seed <- jsonlite::read_json(file.path(outputDirectory, compileInfo))$Rseed[[1]];
                cat("\nRecicling Rseed \"", seed, "\" from last pass\n", sep = "")
            },
            error = function(e) {
                message(e)
                cat(
                    "\nJSON couldn't be read from:\t",
                    compileInfo,
                    ".\t Don't modify this file manually.\n",
                    file = stderr(),
                    sep = "")
            }
        )
    } else {
        seed <- NULL
    }


    layerNames <- preamble_options$layernames
    layerCmd <- preamble_options$layercmd
    answerSheetCorrectTag <- preamble_options$correcttag
    answerSheetWrongTag <- preamble_options$wrongtag

    RandomColumns <- preamble_options$randominfo
    #This shoudl be extra columns with random info.
    # The standard is "randominfo": {"<commandname>": <integer>}}

    ExtraInfo <- preamble_options$extrainfo
    # Extra info will be in the format of "extrainfo" : {"<columnName>":"<cmdName>"}
    #

    sectionReorder <- preamble_options$reordersections

    cmdReorder <- preamble_options$reorderitems


    #### Default behaviour ####


    # Default tables
    if (is.null(tableFile) && is.null(nOutputVersions)) {
        stop("You require a \"table\" of information or a \"noutput\" option to know how many output versions ")
    }

    if (is.null(tableFile)) {
        ClassTable <- data.frame(rollnumber = 1:nOutputVersions, Class = NA, ID = NA, Name = NA, Surname = NA, Nickname = NA)
    } else {
        # We will understand tableFile to be relative to the folder in which the texFile is found, use an absolute path file if you don't want that
        cat(
            tableFile,
            '...exist =',
            assertthat::assert_that(file.exists(tableFile)),
            "\n",
            file = stderr()
        )

        ClassTable <-
            fun_from_folder(
                folder = MainDirectory,
                fun = utils::read.csv,
                tableFile,
                header = TRUE, sep = ",", stringsAsFactors = FALSE,
                na.strings = c("","NaN","NA","Na","NAN"), strip.white = TRUE
            )

        if (is.null(nOutputVersions)) {
            nOutputVersions <- nrow(ClassTable)
        }

    }





    if (is.null(seed)) {
        seed <- as.integer(Sys.time())
    } else if (!is.numeric(seed)) {
        seed <- as.numeric(seed)
        if (is.na(seed)) {
            warning("Seed couldn't be coherced to an integer, using the System time as seed")
            seed <- as.integer(Sys.time())
        }
    }


    if (is.null(layerNames)) {
        layerNames <- c("questions", "choices")

    } else if (!is.character(layerNames)) {
        layerNames <- as.character(layerNames)
    }

    if (is.null(layerCmd)) {
        layerCmd <-  c("question", "(choice|CorrectChoice)")

    } else if (!is.character(layerCmd)) {
        layerCmd <- as.character(layerCmd)
    }

    if (is.null(answerSheetCorrectTag)) {
        answerSheetCorrectTag <- "CorrectChoice"
    } else if (!is.character(answerSheetCorrectTag)) {
        answerSheetCorrectTag <- as.character(answerSheetCorrectTag)
    }

    if (is.null(answerSheetWrongTag)) {
        answerSheetWrongTag <- "choice"
    } else if (!is.character(answerSheetWrongTag)) {
        answerSheetWrongTag <- as.character(answerSheetWrongTag)
    }





    #the seed is set twice, once to generate any random columns
    set.seed(seed)


    # Setting up extra table rows
    #
    #

    #

    cmds <- character()
    columns <- character()

    # Random information
    for (index in seq_along(RandomColumns)) {
        ClassTable[names(RandomColumns)[index]] <-
            sample(
                RandomColumns[[index]],
                size = nrow(ClassTable),
                replace = TRUE
            )

        cmds <- c(cmds, names(RandomColumns)[index])
        columns <- c(columns,  names(RandomColumns)[index])


    }




    ClassTable$rSeed <-
        rep_len(seed, length.out = nrow(ClassTable))


    cmds <- c(cmds, "rseed")
    columns <- c(columns, "rseed")

    # Extra information from the Table
    for (index in seq_along(ExtraInfo)) {

        if (assertthat::is.scalar(ExtraInfo[[index]]) && !is.null(tableFile)) {
            cmds <- c(cmds,  names(ExtraInfo)[index])
            columns <- c(columns, ExtraInfo[[index]])
        }

    }

    #Remove all spaces and punctuations from names, since otherwise it might cause confusion.... keep in mind that R
    #messes with dots and other symbols when reading form a table.

    columns <-
        gsub(
            pattern = "[[:punct:][:space:]]+",
            replacement = "",
            x = columns
        )
    names(ClassTable) <-
        gsub(
            pattern = "[[:punct:][:space:]]+",
            replacement = "",
            x = names(ClassTable)
        )


    # Setting up randomization defaults
    if (!is.logical(sectionReorder)) {
        sectionReorder <- rep_len(TRUE, length(layerNames))
    }
    if (!is.logical(cmdReorder)) {
        cmdReorder <- rep_len(TRUE, length(layerNames))
    } else if (length(cmdReorder) == 1) {
        cmdReorder <- rep(cmdReorder, length(layerNames))
    }

    if (!dir.exists(outputDirectory)) {
        dir.create(outputDirectory)
    }


    #### Calling CreateRandomExams ####


    set.seed(seed) #Reset the seed over here, before randomizing the document

    cat("Generating Exams\n")

    CreateRandomExams(
        x,
        layerNames,
        layerCmd,
        outputBaseName,
        outputDirectory,
        cmdReorder = cmdReorder,
        sectionReorder = sectionReorder,
        infoTable = ClassTable,
        colNames = columns,
        cmdNames = cmds,
        nOutputVersions = nOutputVersions,
        nOutputQuestions = nOutputQuestions,
        answerSheetCorrectTag = answerSheetCorrectTag,
        answerSheetWrongTag = answerSheetWrongTag
    ) ->
        ExamInfo;


    cat("Exams Generated, generating and saving answer sheets...\n", file = stderr())

    ExamSheet <- ExamInfo$FullAnswerSheet

    fileFullAnswerSheet <-
        file.path(
            outputDirectory,
            fullAnswerSheet
        )

    fileShortAnswerSheet <-
        file.path(
            outputDirectory,
            shortAnswerSheet
        )

    #Saving answer sheets for the exam after they have been created
    utils::write.table(
        ExamSheet,
        fileFullAnswerSheet,
        sep = ",",
        quote = FALSE,
        na = "",
        row.names = FALSE
    )


    fileMetaData <-
        file.path(
            outputDirectory,
            compileInfo
        )


    writeLines(
        jsonlite::toJSON(
            c(
                list(Rseed =  seed),
                opt$options,
                preamble_options
                )
        ),
        con = fileMetaData
    )

    OutputShortAnswerTable <- GenerateShortAnswerSheet(ExamSheet)

    ZerothSheet <- stats::na.omit(ExamSheet$CorrectChoice[ExamSheet$Version == 0L])

    if (ncol(OutputShortAnswerTable) == length(ZerothSheet) + 1) {
        # If we can, let's try to add the zeroth version to the answer sheet. (That is, if the number of questions on each exam is the same as the total ammount of questions)

        OutputShortAnswerTable <-
            rbind(
                c(0L, ZerothSheet),
                OutputShortAnswerTable
            )
    }

    utils::write.table(
        OutputShortAnswerTable,
        fileShortAnswerSheet,
        sep = ",",
        quote = FALSE,
        na = "",
        row.names = FALSE
    )

    #### Compile resulting files ####


    if (opt$options$compile) {
        if (opt$options$xelatex) {
            engine <- "xelatex"
        } else {
            # print("Using pdflatex")
            engine <- "pdflatex"
        }
        CompileLatexDir(outputDirectory, outputDirectory, compile.dir = MainDirectory, engine = engine)

        if (!opt$options$debug) {
            #Remove autiliary files

            CompileLatexDir(outputDirectory, outputDirectory, compile.dir = MainDirectory, engine = engine, extracmdoptions = "-c")
        }
    }

}





#' @title Json Homework Parser
#'
#' @description This function takes a series of options as obtained from \code{\link[optparse]{parse_args}} through the parameter \code{opt}. The "examples" section provides all the options that it can parse.
#'
#' From within those options, a \code{--file} option is mandatory.
#'
#' The file option provides a 'LaTeX' file name in which to search for lines on the preamble \code{\%!TexExamRandomizer} within the first 200 lines.
#'
#' With those options that it finds through tags, it passes the function \code{\link{GenerateHomework}}.
#'
#' Note that the tags must respect the JSON format, that is. It \emph{needs} to be written within double quotes.
#' @param opt Options as parsed from \code{\link[optparse]{parse_args}}. The function expects a series of options, the example code exemplifies those options that the function understands.
#'
#' @details
#' It acts similarly to \code{link{jsonexamparser}}, but with the exception of not providing any randomiation option, it only provides the personalization options.
#'
#' Look at \code{vignette("ExamOptions", package = "TexExamRandomizer")} to see the details of the options that it accepts.
#'
#'
#'
#'
#' @family jsoncompiler
#'
#' @export
#' @examples
#' \dontrun{
#' #!/bin/Rscript
#' #This example showcases the type of script this jsonparser might be used on.
#' # You can still use it without a script,
#' # just by adding a list that has the same names as the list provided in opt
#' library(optparse)
#' option_list <- list(
#'     make_option(
#'         c("--file"),
#'         action = "store",
#'         default = NULL,
#'         type = 'character',
#'         help = "Filename of the Tex File"
#'     ),
#'     make_option(
#'         c("--table"),
#'         action = "store",
#'         default = NULL,
#'         type = 'character',
#'         help = "Filename of the table to break down. It overwrites the values written on the file"
#'     ),
#'     make_option(
#'         c("-s", "--seed"),
#'         action = "store",
#'         default = NULL,
#'         type = "integer",
#'         help = "Seed for any randomization done"
#'     ),
#'     make_option(
#'         c("-c", "--compile"),
#'         action = "store_true",
#'         default = FALSE,
#'         type = "logical",
#'         help = "Should the output folder be compiled or not"
#'     ),
#'     make_option(
#'         c("--xelatex"),
#'         action = "store_true",
#'         default = FALSE,
#'         type = "logical",
#'         help = "Should we use xelatex to compile or not"
#'     ),
#'     make_option(
#'         c("-d", "--debug"),
#'         action = "store_true",
#'         default = FALSE,
#'         type = "logical",
#'         help = "If debugging, it doesn't remove auxiliary files"
#'     )
#' )
#'
#'
#' #### PARSING OPTIONS ####
#' ####
#' opt <-
#'     parse_args(
#'         OptionParser(option_list = option_list),
#'         positional_arguments = TRUE
#'     )
#'
#' # Let's assume the file was the example file
#' testfile <-
#'     system.file(
#'         "extdata",
#'         "ExampleTexDocuments",
#'         "exam_testing_nquestions.tex", #Test exam that doesn't require a table
#'         package = "TexExamRandomizer")
#'
#' # To prevent modifying the file system in examples
#' temporalfile <- paste(tempfile(), ".tex", sep = "")
#'
#' file.copy(testfile, temporalfile)
#'
#' opt$options$file <- temporalfile
#'
#'
#'
#' jsonhwparser(opt)
#' }
#'

jsonhwparser <- function(opt) {

    # Checking main tex file exists

    if (is.na(opt$options$file)) {
        stop("Main Tex file not passed")
    } else {
        texFile <- opt$options$file
        cat(
            texFile,
            '...exist =',
            assertthat::assert_that(file.exists(texFile)),
            "\n",
            file = stderr()
        )

    }

    MainDirectory <- dirname(texFile)

    compileInfo <- "compileinfo.txt" # This will hold the information from the compilation steps


    this.basename <- basename(texFile)
    assertthat::assert_that(assertthat::has_extension(texFile, "tex"))
    outputBaseName <- sub(pattern = "\\.tex$", replacement = "", x = this.basename)
    outputDirectory <- file.path(MainDirectory, outputBaseName)


    x <- readLines(texFile)


    # We only search for the %!TexExamRandomizer within the first 200 lines of the document
    preamble_options <- ParsePreambleForOptions(x[1:200])


    if (!is.null(opt$options$table)) {
        tableFile <- opt$options$table
    } else if (!is.null(preamble_options$table)) {
        tableFile <- preamble_options$table
    } else {
        tableFile <- NULL
    }


    if (is.null(opt$options$noutput)) {
        nOutputVersions <- preamble_options$noutput
    } else {
        nOutputVersions <- opt$options$noutput
    }




    if (!is.null(opt$options$seed)) {
        seed <- opt$options$seed
    } else if (!is.null(preamble_options$seed)) {
        seed <- preamble_options$seed
    } else if (file.exists(file.path(outputDirectory, compileInfo))) {
        seed <- NULL

        tryCatch(
            {
                seed <- jsonlite::read_json(file.path(outputDirectory, compileInfo))$Rseed[[1]];
                cat("\nRecicling Rseed \"", seed, "\" from last pass\n", sep = "")
            },
            error = function(e) {
                message(e)
                cat(
                    "\nJSON couldn't be read from:\t",
                    compileInfo,
                    ".\t Don't modify this file manually.\n",
                    file = stderr(),
                    sep = "")
            }
        )
    } else {
        seed <- NULL
    }


    RandomColumns <- preamble_options$randominfo
    #This shoudl be extra columns with random info.
    # The standard is "randominfo": {"<commandname>": <integer>}}

    ExtraInfo <- preamble_options$extrainfo
    # Extra info will be in the format of "extrainfo" : {"<columnName>":"<cmdName>"}
    #



    #### Default behaviour ####


    # Default tables
    if (is.null(tableFile) && is.null(nOutputVersions)) {
        stop("You require a \"table\" of information or a \"noutput\" option to know how many output versions ")
    }

    if (is.null(tableFile)) {
        ClassTable <- data.frame(rollnumber = 1:nOutputVersions, Class = NA, ID = NA, Name = NA, Surname = NA, Nickname = NA)
    } else {
        # We will understand tableFile to be relative to the folder in which the texFile is found, use an absolute path file if you don't want that
        cat(
            tableFile,
            '...exist =',
            assertthat::assert_that(file.exists(tableFile)),
            "\n",
            file = stderr()
        )

        ClassTable <-
            fun_from_folder(
                folder = MainDirectory,
                fun = utils::read.csv,
                tableFile,
                header = TRUE, sep = ",", stringsAsFactors = FALSE,
                na.strings = c("","NaN","NA","Na","NAN"), strip.white = TRUE
            )

        if (is.null(nOutputVersions)) {
            nOutputVersions <- nrow(ClassTable)
        }

    }





    if (is.null(seed)) {
        seed <- as.integer(Sys.time())
    } else if (!is.numeric(seed)) {
        seed <- as.numeric(seed)
        if (is.na(seed)) {
            warning("Seed couldn't be coherced to an integer, using the System time as seed")
            seed <- as.integer(Sys.time())
        }
    }


    #the seed is set twice, once to generate any random columns
    set.seed(seed)


    # Setting up extra table rows


    cmds <- character()
    columns <- character()

    # Random information
    for (index in seq_along(RandomColumns)) {
        ClassTable[names(RandomColumns)[index]] <-
            sample(
                RandomColumns[[index]],
                size = nrow(ClassTable),
                replace = TRUE
            )

        cmds <- c(cmds, names(RandomColumns)[index])
        columns <- c(columns,  names(RandomColumns)[index])


    }




    ClassTable$rSeed <-
        rep_len(seed, length.out = nrow(ClassTable))


    cmds <- c(cmds, "rseed")
    columns <- c(columns, "rseed")

    # Extra information from the Table
    for (index in seq_along(ExtraInfo)) {

        if (assertthat::is.scalar(ExtraInfo[[index]]) && !is.null(tableFile)) {
            cmds <- c(cmds,  names(ExtraInfo)[index])
            columns <- c(columns, ExtraInfo[[index]])
        }

    }

    #Remove all spaces and punctuations from names, since otherwise it might cause confusion.... keep in mind that R
    #messes with dots and other symbols when reading form a table.

    columns <-
        gsub(
            pattern = "[[:punct:][:space:]]+",
            replacement = "",
            x = columns
        )
    names(ClassTable) <-
        gsub(
            pattern = "[[:punct:][:space:]]+",
            replacement = "",
            x = names(ClassTable)
        )


    # Setting up randomization defaults


    set.seed(seed) #Reset the seed over here, before randomizing the document

    #### Generating Homework ####

    if (!dir.exists(outputDirectory)) {
        dir.create(outputDirectory)
    }


    cat("Generating Homework\n")

    GenerateHomework(
        x,
        Table = ClassTable,
        CommandNames = cmds,
        ColumnNames = columns,
        outputDirectory = outputDirectory,
        outputBaseName = outputBaseName
    ) ->
        HWInfo;
    # TODO: Use the filenames from HWInfo to prevent compiling the whole folder, and only compile that that comes frmo HWInfo.
    if (FALSE) { #This is here just to prevent the warning message of "unused variable" until we do the TODO above
        HWInfo
    }



    cat("Homework generated, saving metadata and compiling...\n", file = stderr())


    fileMetaData <-
        file.path(
            outputDirectory,
            compileInfo
        )



    jsonlite::write_json(
        c(
            list(Rseed =  seed),
            opt$options,
            preamble_options
        ),
        path = fileMetaData
    )


    #### Compiling resulting tex files ####

    if (opt$options$compile) {
        if (opt$options$xelatex) {
            engine <- "xelatex"
        } else {
            # print("Using pdflatex")
            engine <- "pdflatex"
        }
        CompileLatexDir(outputDirectory, outputDirectory, compile.dir = MainDirectory, engine = engine)

        if (!opt$options$debug) {
            #Remove autiliary files

            CompileLatexDir(outputDirectory, outputDirectory, compile.dir = MainDirectory, engine = engine, extracmdoptions = "-c")
        }
    }

}
