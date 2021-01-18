# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017

#' @title CountNumberOfSections
#' @description It counts the number of subparts in each section and outputs the result as a table. It doesn't act recursively, it only does the outermost layer.
#' @param Document Document, as defined in \code{\link{StructureDocument}}. Remember however that the function \code{\link{StructureDocument}} returns the document and the preamble together in a list.
#' @details The regular expression that defines this methods behaviour is the following
#'
#' \code{"^[^[:digit:]]*([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)"}
#'
#' The replacement is simply \code{"\\1"}.
#'
#' It tries to first find whether there exist attributes command and section that explain the command and section, before starting to use regexs on the names.
#'
#' @return A table, counting the number of "\code{\\cmdName}" items in which the document was divided when parsed for every begin-end section. It doesn't act resursively.
#'
#' It will return a table with an integer that identifies the section, and a count, with how many items it found on that section. If it doesn't find any items or sections, it will return an empty table.
#' @keywords internal
#' @family Extracting information
CountNumberOfSections <- function(Document){
    #Uses the names instead of greooubg
    #
    sections <- attr(Document, "section")
    commands <- attr(Document, "command")


    if (!is.null(sections)) {
        return(table(sections[sections != 0 & commands != 0 ]))
    }


    #This part down here... preserve until everything is moved to attributes, rather than regexs


    DocNames <- names(Document)
    DocNames <- DocNames[grepl(pattern = "^[^[:digit:]]*([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
                               x = DocNames)]


    sectionNumbers <-
        as.numeric(
            sub(
                pattern = "^[^[:digit:]]*([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
                replacement = "\\1",
                x = DocNames)
        )



    return(
        table(sectionNumbers)
    )
}



#' @title FindExamAnswers
#' @description It finds the answer for a certain document, given a correct and wrong tag. The output character vector is a collection of all those matches, with the text identifying the section and item that it was found inside the tree structure.
#' @details
#' In the document, a correct or wrong item should be identified with a tag. Which shall be a latex command. "\code{\\correct}" "\code{\\wrong}" or whichever. It must be placed somewhere that is not commented-out. (Similar to how the exam class uses the \code{\\CorrectChoice} command to identify a correct answer, instead of using a \code{\\choice}, both tags could be used for this class).
#'
#' (This is internally used by \code{\link{ConstructAnswerSheet}} to construct that DF of answers by then parsing the vector output from this function and getting that way the output)
#'
#' @inheritParams CountNumberOfSections
#' @param correctTag String, it should be the name (or regular expression) defining the tag that items that hold a correct answer will have. \code{\\<correctTag>}.
#' @param wrongTag String, leave as \code{NULL}, unless you went the output to explicitly show those questions that are incorrect. Again, it could be a
#' @param OutputStartingName Internal argument, (it should really be removed after testing that it really does nothing). In theory, this argument starts up the recursive search into the tree for matches. Since the output name will start with this string and a dash afterwards.
#'
#'
#' @return Character vector. Each element identifies one match. The text of the element identifies where that match was found, in terms of the path walked on the tree that it took to get here. The naming convention is specified on details.
#'
#' @details
#' If \code{wrongTag} is not null it also provides the information of where a wrong tag is found.
#'
#' Each output character vectors has a name that identifies all the information necessary to understand where the match was found, relative to both the original document, and the current document we are analizing even if the order of the current document is different.
#'
#' The name of each element starts with \code{<OutputStartingName>}. After that, for each layer that it digs into, it pastes the following name on the right of the name that it already has:
#'
#' \code{_<originalname>_<addedname>}
#'
#' \itemize{
#' \item
#' Where \code{<originalname>} is the name that identifies that list at that depth as the naming convention described in \code{\link{StructureDocument}} defines. It identifies therefore the section and subsection refering to the original document. The string \code{"_original"} is added at the end of the name of the environment and the name of the command name to differentitate it.
#' \item
#' And where \code{<addedname>} is the name that identifies that list at that depth as the naming convention described in \code{\link{StructureDocument}} defines.
#' However, this name has section and item numbers refering to the \strong{current} ordering of the document, not  the original ordering.
#' }
#'
#' In the last layer, when it finally find the correct or wrong tag. It modifies the \code{<addedname>} that should look like \code{i_secName_j_cmdName} and it replaces cmdName with the correct or wrong tag respectively.
#'
#' Therefore, each element is a pretty long string identifying all the layers that it took to traverse to get down to the answer. This function was basically used to prevent the use of attributes that bug out unexpectedly, since when passing functions and parsing things looses the attributes.
#'
#'
#'
#'
#'
#' @keywords internal
#' @family Extracting information

FindExamAnswers <- function(Document, correctTag, wrongTag = NULL, OutputStartingName = "") {

    assertthat::assert_that(is.list(Document))

    matchingVector <-  character()
    pattern1 <- "^[^%]*\\\\"
    pattern2 <- "([^[:alpha:]]|$)"
    patternCorrect <- paste(pattern1, correctTag, pattern2, sep = "")

    patternWrong <-  paste(pattern1, wrongTag, pattern2, sep = "")

    # print(pattern)
    SearchingForTagRecursively <- function(x, level, name.of.x, total.name = name.of.x) {



        if (is.list(x)) {
            q.i <- 0
            s.i <- 0
            current.section.original <- 0
            sections_original <- attr(x, "section_original")
            commands_original <- attr(x, "command_original")
            sections <- attr(x, "section")
            commands <- attr(x, "command")

            for (i in seq_along(x)) {
                if (sections[i] != 0 && commands[i] != 0) {
                    q.i <- q.i + 1
                    # This part over here keeps count of how many questions have we loooked through

                    s.i <- sections[i]

                    # Name relative to the current ordering ->
                    addedName <-
                        sub(
                            "(^[[:digit:]]+)_([^[:digit:]]+?)([[:digit:]]+)_(.*)",
                            paste(s.i, "_\\2", q.i, "_\\4", sep = ""),
                            names(x)[i])

                    # Name relative to the original ordering of the exam, according to the names
                    originalName <-
                        sub(
                            "^([[:digit:]]+_[^[:digit:]]+?)_([[:digit:]]+_[^[:digit:]]+)",
                            "\\1_original_\\2_original",
                            names(x)[i])


                    # Total name to output
                    sub.total.name <-
                        paste(
                            total.name, "_", originalName,"_", addedName, sep = ""
                        )

                    SearchingForTagRecursively(
                        x[[names(x)[i]]],
                        level + 1,
                        names(x)[i],
                        sub.total.name)
                }
            }


        } else if (is.character(x)) {


            if (any(grepl(pattern = patternCorrect, x = x))) {

                matchingVector <<-
                    c(matchingVector,
                      sub("([[:digit:]]+)([^[:digit:]]*?)$",
                          paste("\\1_", correctTag, sep = ""),
                          total.name)
                    )

                # Keep in mind this way of doing it allows more than correctTag match per choice
                # Therefore, maybe when building a table... should keep that in mind
            }else if (!is.null(wrongTag)) { #If wrongTag is null dont try matching thigns
                if (any(grepl(pattern = patternWrong, x = x))) {
                    matchingVector <<-
                        c(matchingVector,
                          sub("([[:digit:]]+)([^[:digit:]]*?)$",
                              paste("\\1_", wrongTag, sep = ""),
                              total.name)
                        )
                }
            }
        }
    }

    SearchingForTagRecursively(Document, 0, OutputStartingName)


    return(matchingVector)
    # This should provide information of the exam, question number, original quesiton number,
    # Choice order, correct choice(Correct choice, only relative to this file)
}


#' @title ConstructAnswerSheet
#' @description Constructs an answer sheet given a document as generated by \code{\link{StructureDocument}} by finding in the items the correct and wrong tags and describing where it found them.
#'
#' Note that you must provide the document part only, \code{StructureDocument} gives back a \code{$preamble} and \code{$document}.
#'
#' If \code{wrongTag} is left \code{NULL}, the answer sheet only shows information of the correct answers.
#'
#' This answer sheet provides information for what answers are correct or incorrect, as well as their position within the original document, before any shuffling was done. (It uses the names of the document to decide whether the document was shuffled or not, since subsetting a list removes all attributes except for the names, this is the "safest" way to do it)
#'
#' The intent of this function is to make it easy to find the answers for a randomized version of an exam.
#'
#' @details
#' The tags are just command of the type "\code{\\Tag}" that must be found somewhere that is not commented out inside the last item at the end of the tree structure. Usually you will want to use the tags that already identify the document items for this.
#'
#' (For example, in the exam class, the tags \code{\\choice} and \code{\\CorrectChoice} could be used naturally, without having to introduce extra commands in the document)
#' @inheritParams CountNumberOfSections
#' @param correctTag Tag to identify the correct items.
#' @param wrongTag Tag that identifies the wrong items.
#' @return Data Frame. With the following columns
#' \describe{
#'     \item{index}{Just an index running from 1 to \eqn{n}, where \eqn{n} is the numbe of rows}
#'     \item{For each layer of depth in the document:}{
#'     Four columns,
#'     \describe{
#'         \item{\code{<name of section>_original}}{Contains an integer identifying the numbering of this section in the original layer, as identified by the naming convention}
#'         \item{\code{<name of section command>_original}}{Contains an integer identifying the numbering of this item in the original section, as identified by the naming convention}
#'         \item{\code{<name of section>}}{Contains an integer identifying the numbering of this section in the current layer, as identified by the ordering of the  document inputted on this function}
#'         \item{\code{<name of section command>}}{Contains an integer identifying the numbering of this item in the current section, as identified by the ordering of the  document inputted on this function}
#'     }
#'     }
#'     \item{For the last layer of depth}{
#'     5 columns if the wrongTag is not NULL, 4 columns otherwise,
#'     \describe{
#'         \item{\code{<name of section>_original}}{Contains an integer identifying the numbering of this section in the original layer, as identified by the naming convention}
#'         \item{\code{<name of section command>_original}}{Contains an integer identifying the numbering of this item in the original section, as identified by the naming convention}
#'         \item{\code{<name of section>}}{Contains an integer identifying the numbering of this section in the current layer, as identified by the ordering of the  document inputted on this function}
#'         \item{\code{<correctTag>}}{
#'             Contains an integer identifying the numbering of this item in the current section, , as identified by the ordering of the  document inputted on this function
#'
#'             If the \code{correctTag} wasn't found in this item, it will show \code{NA} instead. (This will only happen if \code{wrongTag} is not \code{NULL}, since otherwise this elements are omitted)}
#'         \item{\code{<wrongTag>}}{
#'             Contains an integer identifying the numbering of this item in the current section, as identified by the ordering of the  document inputted on this function
#'
#'             If the \code{wrongTag} wasn't found in this item, it will show \code{NA} instead. (This will only happen if \code{wrongTag} is not \code{NULL}, since otherwise this elements are omitted)}
#'     }
#'     }
#' }
#'
#' @seealso \code{\link{FindExamAnswers}} for the exact underlying messy algorithm that controls how the table is created.
#' @family Extracting information
#' @export
#' @examples
#'
#' ConstructAnswerSheet(
#'     TexExamRandomizer::testdoc$document,
#'     "CorrectChoice",
#'     "choice"
#' )

ConstructAnswerSheet <- function(Document,
                                 correctTag,
                                 wrongTag = NULL) {




    ExamAnswers.char.vector <- FindExamAnswers(Document, correctTag, wrongTag)
    # print(ExamAnswers.char.vector)
    result <- data.frame(index = seq_along(ExamAnswers.char.vector))

    while (any(grepl("^[^[:digit:]]*[[:digit:]]+_[^[:digit:]]+?_[[:digit:]]+_[^[:digit:]]*?_",
                     ExamAnswers.char.vector)
    )) {
        #Outer Layers
        sectionNumberOuter <-
            as.integer(
                sub("^[^[:digit:]]*([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
                    "\\1",
                    ExamAnswers.char.vector
                )
            )

        TagNumberOuter <-
            as.integer(
                sub(
                    "^[^[:digit:]]*([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
                    "\\3",
                    ExamAnswers.char.vector
                )
            )
        sectionOuter <-
            sub(
                "^[^[:digit:]]*([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
                "\\2",
                ExamAnswers.char.vector
            )
        TagOuter <-
            sub(
                "^[^[:digit:]]*([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_([^[:digit:]]*)_(.*)",
                "\\4",
                ExamAnswers.char.vector
            )


        ExamAnswers.char.vector <-
            sub(
                "^[^[:digit:]]*[[:digit:]]+_[^[:digit:]]+?_[[:digit:]]+_[^[:digit:]]*?_",
                "",
                ExamAnswers.char.vector
            )


        result <- cbind(result, sectionNumberOuter, TagNumberOuter, stringsAsFactors = FALSE)

        n <- length(result)
        names(result)[(n - 1):n] <- c(sectionOuter[1], TagOuter[1])



        #Make sure to go back and disallow numbers on names of thigns back in the structure functions to make sure this funciton doesn't break
    }

    #Answer's Layer
    # The last layer has a different structure
    sectionNumberInner <-
        as.integer(
            sub(
                "[^[:digit:]]*?([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
                "\\1",
                ExamAnswers.char.vector
            )
        )
    TagNumberInner     <-
        as.integer(
            sub(
                "[^[:digit:]]*?([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
                "\\3",
                ExamAnswers.char.vector
            )
        )
    sectionInner <-
        sub(
            "[^[:digit:]]*?([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
            "\\2",
            ExamAnswers.char.vector
        )
    TagInner <-
        sub(
            "[^[:digit:]]*?([[:digit:]]+)_([^[:digit:]]+)_([[:digit:]]+)_(.*)",
            "\\4",
            ExamAnswers.char.vector
        )

    # LastTag, WHY:
    # There could be more than one possible final tag, since there is a correct and wrong tag (And maybe more if we extend the funcitonality to find extra tags)
    # We want to divide this last column into different columns. Each column will identify the numbers for only the tag that represents them and will have NA on all the different entries.
    # For example, one column could be the correctChoice column, and identifies the number of the last choice, if the naswer is correct or NA if that row is not identifying a correct choice. And the next column couldbe the wrong choices.
    allTags <- stats::na.omit(unique(TagInner))
    result <- cbind(result, sectionNumberInner, stringsAsFactors = FALSE)

    for (Tag in allTags) {
        extracolumn <- TagNumberInner
        extracolumn[TagInner != Tag] <- NA


        result <- cbind(result, extracolumn, stringsAsFactors = FALSE)

    }


    n <- length(result)
    numTags <- length(allTags)

    names(result)[(n - numTags):n] <- c(sectionInner[1], allTags)

    return(result)


}


#### TODO5?????: Allow to have choices removed as well? but always keeping the correctchoice in there
#### huh? I wrote that up there way too long ago, can't rememeber what TODO5 was meant to say
# answer.DF <- ConstructAnswerSheet(b, "CorrectChoice", "choice")
#
#



#' Generating a short answer sheet
#'
#' Given a number of answer sheets generated by \code{\link{ConstructAnswerSheet}} that have been binded together. And that have a column, \code{versionColName}, that identifies each version. It collects all the answers together and places all the answers together for each exam.
#'
#' Note that if the version number is 0, it is ignored, since it understands that version 0 is the reference version.
#'
#' If the document has more than two layers, keep in mind that it just shows the top most layer numbering and then the inner most number of the correct answers.
#'
#' Note how this implies as well that an exam with more than one possible answer can not be simplified into a short answer sheet.
#'
#' @param ExamSheet a exam sheet that contains all versions, similar to \code{\link{CreateRandomExams}}
#' @param versionColName The name of the column in the original exam that contains the version number
#' @param correctColName The name of the column that contains the last index for the correct tag, or NA if it is not a correct choice.
#'
#'
#' @return A data frame \itemize{
#'  \item Each row identifies one version of the answer sheet
#'  \item the first column is the version number, the rest of the columns are the questions,
#'  }
#' @details
#' IMPORTANTLY, If a certain exam has less answers than other exams, the are just cited sequentially. Which may cause confusion. To clarify. This may happen if a certain question has more than one solution marked as "correct", or if a certain question has no solutions marked as correct. In that case, The short answer sheet just sequentially names all the correct answers, disregarding which questions they are referring to. (This is a very special case that will only come up in a real scenario if you are writing a short answer question in the middle of a multiple choice test. Or if you are writing some questions to have multiple correct answers, but only a few of them, and those questions are not included in all exams... (So evil))
#' @export
#' @family Extracting information
#' @examples
#'
#' csvfile <- system.file(
#'     "extdata",
#'     "ExampleTables",
#'     "ExampleAnswerSheet.csv",
#'     package = "TexExamRandomizer"
#' )
#' testASheet <- read.csv(
#'     csvfile,
#'     header = TRUE,
#'     stringsAsFactors = FALSE,
#'     na.strings = c("", "NA", "Na"),
#'     strip.white = TRUE
#' )
#'
#' GenerateShortAnswerSheet(testASheet)

GenerateShortAnswerSheet <- function(ExamSheet, versionColName = "Version", correctColName = "CorrectChoice") {

    DF <-
        ExamSheet[
            ExamSheet[[versionColName]] != 0L,
            ]
    VersionList <- unique(DF[[versionColName]])
    v.1stLine <- #Vector of all the correct choices without NA,s
        stats::na.omit(
            DF[
                DF[[versionColName]] == VersionList[1L],
                ][[
                    correctColName
                    ]]
        )
    nQuestions <- length(v.1stLine) # Since that is the correct
    nVersions <- length(VersionList)
    result <- matrix(ncol = nQuestions + 1L, nrow = length(VersionList))
    result[1,] <-
        c(
            VersionList[1L],
            v.1stLine
        )
    if (nVersions > 1) {
        for (i in 2:nVersions) {
            answer_row <-
                stats::na.omit(
                    DF[
                        DF[[versionColName]] == VersionList[i],
                        correctColName
                        ]
                )
            if (length(answer_row) < nQuestions) {

                answer_row <- c(answer_row, rep(NA, nQuestions - length(answer_row)))

            } else if (length(answer_row) > nQuestions) {

                result <- cbind(result, matrix(data = NA, nrow = nrow(result), ncol = length(answer_row) - nQuestions))
                nQuestions <- length(answer_row)

            }


            result[i, ] <-
                c(
                    VersionList[i],
                    answer_row
                );
        }
    }
    return(
        structure(
            as.data.frame(result),
            names = c(
                versionColName,
                paste(correctColName, 1:nQuestions)
            )
        )
    );
}

