# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017



#### FUNCTIONS ####
#' @title Generate Homework
#' @description This function personalizes a 'LaTeX' document with data from a table,
#' generating a new file for each row which is saved on the \code{outputDirectory}.
#'
#' @details The command names should be 'LaTeX' commands that are being defined through
#'
#'  \code{\\newcommand\{\\<CommandNames[i]>\}\{<previous definition>\}}
#'
#'  The definition of these commands will be changed to be
#'
#'  \code{\\newcommand\{\\<CommandNames[i]>\}\{<Table[ColumnNames[i]][file #]>\}}
#'
#'  And it will output one file for each command.
#'
#'  The intent of this function was to populate information into a generic homework to personalize it for every student using 'LaTeX'.
#'  (It actually generalizes to maybe other problems).
#' @inheritParams DivideFile
#' @param Table Data frame from which to extract the information
#' @param CommandNames Character vector with the same length as \code{columnNames}
#' @param ColumnNames Character vector with the names of the columns to be used
#' @param outputDirectory The directory in which the output will be placed
#' @param outputBaseName The starting name for the output files
#'
#' The files will look like
#'
#' \code{<outputDirectory>/<outputBaseName>_00<number>.tex}
#'
#' Where the number of zeros is the minimum number of zeros required to have a different version number for each file. (i.e., if there is only 45 files, it is 01-45; but with 132 files, it would be 001-132)
#'
#' @seealso \code{\link{ReplaceFromTable}} to get a better idea of how the replacement is made. To see examples of how to use it, look at the code in \code{\link{jsonhwparser}}
#' @return Character vector with the file names of the output.
#'
#'
GenerateHomework <- function(
  x,
  Table,
  CommandNames,
  ColumnNames,
  outputDirectory,
  outputBaseName
) {


  #Negative numbers will work incorrectly
  # If outputcharlen < length(int.vector, it will throw an error)
  # Returns a character vector adding to the int.vector 0 until it has the same character length as outputcharlen asks for.
  #

  NumberwithZeros <- function(int.vector, outputcharlen){
    result <- character(length = length(int.vector))
    for (k in seq_along(result)) {

      result[k] <-
        paste(
          paste(
            rep("0", outputcharlen - nchar(int.vector[k])),
            collapse = ""
          ),
          int.vector[k],
          sep = ""
        )

   }
    return(result)
  }

  # Outputs using cat a tex document given it's preamble, maindocument and the filename to output it to.
  # NOTE: After each character vector element is outputted to the file a \n is added.
  # Args:
  #  preamble: first list or character vector. It will contain the first part of the document.
  #  document: second list or character vector, will be outputted to the stream after preamble
  #  this.filename: filename in which to stream, a .tex is added in the end
  #
  # Returns:
  #  the total file path in which the stream was outputted in which it was outputted
  SaveHomework <- function(preamble, document, this.filename){


    file <- paste(this.filename, ".tex", sep = "")

    catDocument(list(preamble, document), file = file)

    return(file)

  }



  Document <- DivideFile(x)

  nOutputFiles <- nrow(Table)

  outputFileNames <-
    paste(
      file.path(outputDirectory,
                outputBaseName
      ),
      "_",
      NumberwithZeros(
        1:nOutputFiles,
        nchar(nOutputFiles)
      ),
      sep = ""
    );
  # The .tex is added by SaveHomework... don't ask me why, it seems silly, but if later on some other formatting was needed, it would be easy to just change the SaveHomework function instead of digging here


  #IF not all the columnnames are in the table, try to use some regex logic to match them:
  if (!all(ColumnNames %in% names(Table))) {
    TableNames <- names(Table)
    for (i in seq_along(ColumnNames)) {
      if (!(ColumnNames[i] %in% TableNames)) {
        tryCatch(
          {
            match <-
              match(
                TRUE,
                grepl(
                  pattern = ColumnNames[i],
                  x = TableNames,
                  ignore.case = TRUE
                )
              )
          },
          error = function(e) stop("Column name not found and can't be used in regex")
        )
        # The main error we are trying to catch here is if the regex pattern of colname is "incorrect" (syntactically).
        # Function should stop in that case

        if (assertthat::assert_that(!is.na(match))) {
          # If there was a match
          ColumnNames[i] <- TableNames[match]
        }
      }
    }
  }
  #assert_that(all(ColumnNames %in% names(Table)))


  for (i in 1:nrow(Table)) {

#    cat("Generating Homework",i, "\nInfo:", Table[ColumnNames][i, ] %>% unlist, "\n")

    cat(
      "Generating Homework",i,
        # "\nInfo:", unlist(Table[ColumnNames][i, ]),
        "\n",
      file = stdout() # should this be outputted on the stderr, which could be captured in a log. ??
      )

    OutputPreamble <-
      ReplaceFromTable(
        x = Document[[1]],
        table = Table,
        tableRow = i,
        columnNames = ColumnNames,
        commandNames = CommandNames
      );

    cat("\tFinished with Homework",i, "\n\t\t... Saving it to:\n\t\t...", outputFileNames[i], "\n\n")


    SaveHomework(
      OutputPreamble,
      Document[[2]],
      this.filename = outputFileNames[i]
    )



  }


  return(outputFileNames)

}




#' @title CreateRandomExams
#'
#' @description
#' This function creates a series of randomized exams from a tex document and personalizes the information from a table (if a table is given) and a series of command names where thae information shoudl be replaced.
#'
#' @details
#' All the output exams are named with \code{outputBaseName} followed by 00i identifying the number of the exam (The number of zeros is the minimum that allows for all the exams to have a different number) and \code{"_Version_"} followed by the version number of the exam and "\code{.tex}". That is:
#'
#' \code{<outputDirectory>/<outputBaseName>00i_Version_j.tex}
#'
#' The number of exams outputted will always be the same as the number of versions if no table is given. However, if a table is added as input. It will create one exam for each row of the table, and it will try to divide as evenly as possible how to give the versions between the different rows. (Having one exam for each row, which will probably represent a student)
#'
#'
#' @inheritParams StructureDocument
#' @param outputBaseName String, The basename for the output files.
#' @param outputDirectory String, The output directory.
#' @param cmdReorder,sectionReorder Logical vector, the length of \code{cmdReorder} determines how many layers deep are we going to dig and randomize. For that reason, if \code{sectionReorder} is just a scalar, it will assume that it repeats for every \code{cmdReorder} that is given. See \code{\link{RandomizeDocument}} for extra details on these parameters.
#' @param infoTable Table with information, if NULL, no information is added to the exams
#' @param colNames Character vector, Column names from the \code{infoTable} from which we will extract the information.
#'
#'        It first tries to find the column names literally, if ti couldn't find them like that, it will try to use them as a regular expression to find a column that matches the column.
#' @param cmdNames Character vector, Names of the commands on the tex file, \code{\\<cmdNames[i]>}, that are to be matched with the columns to replace the information from the table in those commands. For extra info see also \code{\link{ReplaceFromTable}}
#' @param nOutputVersions Number of different random versions of the exam to be outputted
#' @param nOutputQuestions Number of "questions" on the output exams. If the input is a scalar, the program will decide how to more evenly split the questions between all the sections, otherwise one can directly provide an integer vector specifying how many questions from each section are needed. (this only searches the "items" of the outermost layer)
#' @param answerSheetCorrectTag,answerSheetWrongTag If the tags are not given, the output answersheet will be \code{NULL}. In other cases, these tags can be regular expressions
#' @param optionList Instead of writing the options on the function. Options could be given to optionList, and it will add those options. As long as the names are correct
#' @return A list that contains
#' \describe{
#'     \item{\code{outputDirectory}}{The output directory}
#'     \item{\code{outputFiles}}{A character vector that contains all the output names}
#'     \item{\code{FullAnswerSheet}}{The full answer sheet of all the exams.
#'
#'     Each answer sheet is created as described by \code{\link{ConstructAnswerSheet}}, and all the answer sheets are joined together with a version number in front as an added column to bind them all together. The original version has the number 0, all the output versions have sequential numbers as Version
#'
#'     This wrapper function assumes equal depth on all branches of the tree structure, so that the number of columns is always identical in the answer sheet
#'     }
#' }
#'
#' @seealso \code{\link{ConstructAnswerSheet}}, \code{\link{ReplaceFromTable}}, \code{\link{RandomizeDocument}} for extra details. . To see examples of how to use it, look at the code in \code{\link{jsonhwparser}}
#' @export

CreateRandomExams <- function(
    x,
    layersNames = c("questions", "choices"), layersCmd = c("question", "(choice|CorrectChoice)"),
    outputBaseName,
    outputDirectory,
    cmdReorder = rep_len(TRUE, length(layersNames)),
    sectionReorder = FALSE,
    infoTable = NULL,
    colNames = NULL,
    cmdNames = NULL,
    nOutputVersions = nrow(infoTable),
    nOutputQuestions = "max",
    answerSheetCorrectTag = NULL,
    answerSheetWrongTag = NULL,
    optionList = NULL
){

    if (!is.null(optionList)) {
        for (option_number in seq_along(optionList)) {
            assign(names(optionList)[option_number], optionList[[option_number]])
        }
    }


    # This function si a wrapper for creating from a tex file a bunch of randomly generated exams, following the structure of a list of lists as described in the documentation of this packet.
    #
    # This function will save a bunch of tex files in the output directory given. All of them being reordered and with personalized information versions of the starting exam given.
    #
    # Args:
    # x: The file as a character vector (One line of file per line of the vector)
    # layersNames: Names of the layers sectioning \begin{layersName} etc
    # layersCmd: Names of the "item" commands \layersCmd
    # outputBaseName: Name of the base that will be used to generate exams.
    # outputDirectory: Directory in which the output will be stored
    # cmdReorder: Should the layer be reordered or not,
    # sectionReorder: Should the different \begin \end sections reordered or not
    # infoTable: Table from which to input extra information into the preamble, (Will be saved in commands)
    # colNames: Names of the columns form the table that contain the information. Regex expressions can be used for a more flexible naming convention on the csv tables
    # cmdNames: Names of the commands on the latex file that those columns contain information to (Basically it searches for \newcommand{\cmdName}{bla bla} and replaces the bla bla for any information)
    # nOutputVersions: Number of output versions of the exam, the predetermined value is one for each student on the table if a table is given
    # nOutputQuestions: How many questions should it be. If it can't be casted as integer, or NA, or NULL, it will ouput an exam with the maximum number of questions possible.
    # answerSheetCorrectTag: If generating and answer sheet, what is the cmd that marks an answer to be correct. Generating an answer sheet is activated when the value of this variable is not null
    # answerSheetWrongTag: If generating an answer sheet, what is the cmd of an answer that is wrong
    # Returns:
    # list...
    #   ...$outputDirectory: outputDirectory,
    #   ...$outputFiles: outputFileNames,
    #   ...$FullAnswerSheet: The full answer sheet for all the exams
    #     As described on the function ConstructAnswerSheet,
    #     the structure is fairly flexible.
    #

    if (length(sectionReorder) == 1) {
        sectionReorder <- rep_len(sectionReorder, length(cmdReorder))
    }

    #Determining if we are creating an answer sheet to be returned.
    if (is.null(answerSheetCorrectTag)) {
        createAnswerSheet <- FALSE
    } else {
        createAnswerSheet <- TRUE
    }


    # This functions returns the indices to be removed from a sectionlist to return the indices needed to be removed.
    # Avoids removing questoins
    #
    # Then one can just do a - in front of it to remove them
    #
    # Args:
    #   SectionNames: Character vector, Names of the list to be understood,  following the convention of the "texTreefication" functions
    #   remove.vector: An integer vector with as many items as sections indicating how many quesitons should be removed, it removes things from the end back.
    #
    # Returns:
    #   An integer vector with the indices of the items to be removed. Do a list[-result] to keep those that you want.
    RemoveLastQuestions <- function(SectionNames, remove.vector){

        assertthat::assert_that(is.character(SectionNames))

        # result.which announces which indices to remove, we are doing the opposite, since it is easier to remove an amount of elements m, instead of adding n elements.

        result.which <- integer()

        for (i in seq_along(remove.vector)) {
            result.which <-
                c(
                    utils::tail(
                        which(
                            grepl(
                                pattern = paste(
                                    "^[^[:digit:]]*",
                                    i,
                                    "_([^[:digit:]]+)_([[:digit:]]+)",
                                    sep = ""
                                ),
                                x = SectionNames
                            )
                        ),
                        remove.vector[i]
                    ),
                    result.which
                )

        }


        return(result.which)


    }

    # This function tales the nOutputQuestions that we want and returns the remove vector, detailing for each section how many parts should be removed to match the nOutputquestions requested, removing the questions as evenly as possible
    #nOutputQuestiosn should be an integer vector
    #FullDocumentNumberQ should be a vector saying for each section of the main document, how many questions there are
    #
    # This function returns 0 if there is no need to remove vector

    GetRemoveVectorFromInput <- function(nOutputQuestions, FullDocumentNumberQ){

        assertthat::assert_that(is.numeric(nOutputQuestions), is.numeric(FullDocumentNumberQ))

        if (length(nOutputQuestions) != 1L && length(nOutputQuestions) != length(FullDocumentNumberQ)) {
            stop("Incorrect number of sections detected. The \"number of questions\" vector should have the same length as the number of sections on this document, or be a scalar and allow the program to distribute where to remove the questions between the sections.")
        }


        total.OutputQuestions <- sum(nOutputQuestions) # Exam output
        total.MainQuestions   <- sum(FullDocumentNumberQ) #Input Exam total number
        nSections <- length(FullDocumentNumberQ)

        # We need to differentiate between when the nOutput is given as one element or as a vector

        # Set up for the function
        if (length(nOutputQuestions) == 1) {
            if (total.OutputQuestions < nSections) {

                warning("You want an exam with less questions than sections, palomar principle tells me that I can't, the minimum is one question per section. Which is what I will do")

                removeQuestions.vector <- pmax.int(FullDocumentNumberQ - 1, 1)
            } else if (total.OutputQuestions >= total.MainQuestions) {
                warning("Output number of questions that was requested is greater than the number of questions we can count on the main file")
                return(0L)

            }


        } else {# in this case nOutput is a vector
            if (any(nOutputQuestions < 1)) {

                nOutputQuestions <- pmax.int(nOutputQuestions, 1)
                warning("Each section must have at least one part. Otherwise, just remove the section all together from the original file, setting the output to be one question per section, if this is not the desired result, change the settings.")

            }

            if (any(nOutputQuestions > FullDocumentNumberQ)) {
                warning("You are giving a number of  question that was requested for certain sections greater than those we can count on the file")
                nOutputQuestions <- pmin.int(nOutputQuestions, FullDocumentNumberQ)
            }

        }

        # Function itself. we have two cases, if the input is just one element, or if it is a vector
        #


        if (length(nOutputQuestions) == 1) {

            removeQuestions <- total.MainQuestions - total.OutputQuestions

            removeQuestions.vector <- #This vector says for each section how many quesiotns should be removed
                floor(
                    rep(
                        removeQuestions, nSections
                    ) * FullDocumentNumberQ / total.MainQuestions
                )


            # And we may need to get those that we missed

            while (sum(removeQuestions.vector) < removeQuestions) {

                #The rest of the questions, we remove them from random sections
                # Trying to guarantee we always leave one, and sections that have more questions are more likely to get something removed.
                # #Note how all exams will have the same number of questions removed for each section, but not necesarily the same questions will be removed in each section, that part will be random.

                # Getting the random index frmo which to remove
                k <-
                    sample.int(
                        nSections,
                        size = 1,
                        prob = pmax(FullDocumentNumberQ - removeQuestions.vector - 1, 0)
                    );

                removeQuestions.vector[k] <- removeQuestions.vector[k] + 1
            }

            # Removes from each section proportionally to the amount of questions in that section, all exams will have the same amount of questions per section.
            #

        } else {

            removeQuestions.vector <- FullDocumentNumberQ - nOutputQuestions

        }

        return(removeQuestions.vector)

    }

    # This function takes integers and returns them written as "000int" with as many zeros as needed to make all the vectors the same length character wise
    #Negative numbers will work incorrectly
    # Args:
    #   int.vector: Integer vector with all the integers to be transformed
    #   outputcharlen: Output digit length that we require, should be bigger than the biggest number
    # If outputcharlen < length(int.vector, it will throw an error)
    #
    # Returns:
    #   Character vector with "000<int>" for each int in the input vector, where the number of characters of each outputcharacter is outputcharlen
    NumberwithZeros <- function(int.vector, outputcharlen){


        result <- character(length = length(int.vector))
        for (k in seq_along(result)) {
            result[k] <-
                paste(
                    paste(
                        rep("0", outputcharlen - nchar(int.vector[k])),
                        collapse = ""
                    ),
                    int.vector[k],
                    sep = ""
                )
        }
        return(result)
    }


    # Saves a tex file into the file given, following the patern
    # <this.filename>_Version_<VersionNumber>.tex
    # Uses the function cat to do so.
    #
    # Args:
    #   preamble: Preamble of the tex file (from start to \begin{document})
    #   document: Document of the tex file (from \begin{document} until \end{document})
    #   this.filename: desired output file name
    #   versionNumber: Output version
    #
    # Returns:
    #   Returns the address of the file in which the file was saved

    # Returns the file in which it was outputted

    SaveExam <- function(preamble, document, this.filename, versionNumber){
        file <- paste(this.filename, "_Version_", versionNumber, ".tex", sep = "")

        cat("Saving Exam with version", versionNumber, "...")

        catDocument(list(preamble, document), file = file)

        cat("\tDone\n")


        return(file)

    }







    TexDocument <-
        StructureDocument(
            x,
            layersNames = layersNames,
            layersCmd = layersCmd
        )

    # Determining if we need to remove question to achieve the desired result
    #
    # Checking for nOutputQuestions
    nOutputQuestions <- as.integer(nOutputQuestions)
    if (is.na(nOutputQuestions) || length(nOutputQuestions) < 1) {
        # If we don't understand the nOutputQuestions properly, don't try to remove questions, assume the number of questions is the maximum
        RemoveQuestions <- FALSE
    } else {
        FullDocumentNumberQ <-
            as.numeric(
                CountNumberOfSections(TexDocument[[2]])
            )

        removeQuestions.vector <- GetRemoveVectorFromInput(nOutputQuestions, FullDocumentNumberQ)

        RemoveQuestions <- any(removeQuestions.vector > 0L)

    }

    if (createAnswerSheet) {
        #Initialize the output answer sheet with the barebones input exam.
        answer.DF <-
            cbind(
                Version = 0,
                ConstructAnswerSheet(
                    TexDocument[[2]],
                    answerSheetCorrectTag,
                    answerSheetWrongTag
                )
            )

    } else {
        # Otherwise it will output NULL
        answer.DF <- NULL
    }

    #Checking for nOutputVersions
    #
    nOutputVersions <- as.integer(nOutputVersions)
    if (is.null(nOutputVersions) || any(is.na(nOutputVersions)) || length(nOutputVersions) != 1) {
        #If we don't know what it is
        if (!is.null(infoTable)) {
            nOutputVersions <- nrow(infoTable)
        } else {
            stop("The number of output versions wasn't specified properly")
        }
    }

    # Checking for info Table
    if (is.null(infoTable)) {
        nOutputFiles <- nOutputVersions

    } else {

        nOutputFiles <- nrow(infoTable) # One line per file

        if (nrow(infoTable) > nOutputVersions) {
            # If we have more rows than versions, we have to repeat versions for each student

            nFilesPerVersion <-
                rep(nOutputFiles, nOutputVersions) %/% nOutputVersions
            # Divide integer wise to somewhat evenly distribute it
            #


            # TODO: Improve the way we are calculating this... this part is kind of a mess


            # Add the rest of exams that we need left
            # nFilesPerVersion[i] will be the number of different outputs we will have for each exam with version [i]
            j <- 1
            while (sum(nFilesPerVersion) < nrow(infoTable)) {
                # If we can't evenly distribute the exam versions between the different students on the table, extra exams are given in sequential order.
                nFilesPerVersion[j] <- nFilesPerVersion[j] + 1
                #cycle around
                if (j == length(nFilesPerVersion)) {
                    j <- 1
                } else {
                    j <- j + 1
                }
            }

        } else {
            nFilesPerVersion <- rep(1, nOutputVersions)
        }

        ## And now that everything is set up let's do the regex to find the columsn that are needed.
        ##

        if (!all(colNames %in% names(infoTable))) { #Try to regex the names from the columns into the actual table if they were not all found exactly
            TableNames <- names(infoTable)
            for (i in seq_along(colNames)) {
                if (!(colNames[i] %in% TableNames)) {
                    tryCatch(
                        {
                            match <-
                                match(
                                    TRUE,
                                    grepl(
                                        pattern = colNames[i],
                                        x = TableNames,
                                        ignore.case = TRUE
                                    )
                                )
                        },
                        error = function(e) stop("Column name not found and can't be used in regex")
                    )

                    assertthat::assert_that(
                        !is.na(match),
                        msg = paste("Column pattern \"", colNames[i], "\" was not found. Haulting execution", sep = "")
                        )


                    #If there was a match
                    colNames[i] <- TableNames[match]
                }
            }
        }

    }





    outputFileNames <-
        paste(
            file.path(outputDirectory,
                      outputBaseName
            ),
            NumberwithZeros(
                1:nOutputFiles,
                nchar(nOutputFiles)
            ),
            sep = ""
        ); #the names of the output files without the Version # and the .tex, that will be added later

    currentExamNumber <- 0

    for (version.i in seq_len(nOutputVersions)) {

        cat("Generating Exam ", version.i, ":...\n", sep = "")

        cat("\tAdding version...")

        version.preamble <-
            ReplacePreambleCommand(
                x = TexDocument[[1]],
                commandName = "myversion",
                commandValue = version.i
            )
        cat("\tDone\n")


        cat("\tRandomizing Exam...")

        rndDocument <-
            RandomizeDocument(
                Document = TexDocument[[2]],
                isSectionReordered.vector = sectionReorder,
                isLayerRandomized.vector = cmdReorder
            )
        cat("\tDone\n")

        if (RemoveQuestions) {
            remove_questions <- RemoveLastQuestions(names(rndDocument), removeQuestions.vector)

            rndDocument <- SubsetWithAtributes(rndDocument, -remove_questions)
        }

        if (createAnswerSheet) {

            answer.DF <-
                rbind(
                    answer.DF,
                    cbind(
                        Version = version.i,
                        ConstructAnswerSheet(
                            rndDocument,
                            answerSheetCorrectTag,
                            answerSheetWrongTag
                        )
                    )
                );

        }


        if (!is.null(infoTable)) {
            for (j in seq_len(nFilesPerVersion[version.i])) {
                currentExamNumber <- currentExamNumber + 1
                preamble <-
                    ReplaceFromTable(
                        x = version.preamble,
                        table = infoTable,
                        tableRow = currentExamNumber,
                        columnNames = colNames,
                        commandNames = cmdNames
                    )
                # ReplaceFromTable <- function(x, table, tableRow, columnNames, commandNames)
                #

                outputFileNames[currentExamNumber] <-
                    SaveExam(
                        preamble = preamble,
                        document = rndDocument,
                        this.filename = outputFileNames[currentExamNumber],
                        versionNumber = version.i
                    )
            }

        } else {

            currentExamNumber <- currentExamNumber + 1
            # ReplacePreambleCommand <- function(x, commandName, commandValue)
            #
            outputFileNames[currentExamNumber] <-
                SaveExam(
                    preamble = version.preamble,
                    document = rndDocument,
                    this.filename = outputFileNames[currentExamNumber],
                    versionNumber = version.i
                )

        }
    }



    return(
        list(
            outputDirectory = outputDirectory,
            outputFiles = outputFileNames,
            FullAnswerSheet = answer.DF
        )
    );


}


