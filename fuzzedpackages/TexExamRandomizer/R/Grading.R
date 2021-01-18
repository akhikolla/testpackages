# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017
#### SET UP FUNCTIONS ####
####
ORIGINAL_PATTERN <- "_original$"
#just to remember the convention of how an original column looks on the original file


#' @title FindMatchingRow
#'
#' @description
#' It outputs a logical vector identifying which rows in \code{DF} have the same values as the values on \code{rowtoMatch}. It ignores all columns that are not found on \code{rowtoMatch} when doing the matching.
#'
#' @param rowtoMatch One row of a data frame. If the length is longer, it won't output it.
#' @param DF Dataframe in which we want to find the matches.
#'
#' It doesn't search that both the row of Df and rowtoMatch matches exactly, it only checks whether the columns that are found on rowtoMatch are found on DF with the same values (It decides which columns are "the same column" by looking at the names of the columns).
#'
#' It will throw an error if the names of the columns on \code{rowtoMatch} are not all found on \code{DF}.
#' @return Logical vector, with the same length as the number of rows in \code{DF}. It outputs \code{TRUE} if that row matched the \code{rowtoMatch}, \code{FALSE} otherwise
#'
#' @details
#' \strong{Note:} It matches using the name of columns from \code{rowtoMatch}, matching them to the names of the \code{DF}.
#'
#'
#' It will output \code{logical(0)} if the \code{rowtoMatch} is null, or if is a data frame with 0 rows.
#'
#'
#' @keywords internal
#' @family Grading exams core
FindMatchingRow <- function(rowtoMatch, DF) {

  if (is.null(rowtoMatch)) {
    warning("NULL rowtoMatch")
    return(logical(0))
  } else if (nrow(rowtoMatch) == 0) {
    warning("no rows to match, 0 rows")
    return(logical(0))
  }
  assertthat::assert_that(all(names(rowtoMatch) %in% names(DF)))
  assertthat::assert_that(nrow(rowtoMatch) == 1)

  result <- rep(TRUE, nrow(DF))

  for (colName in names(rowtoMatch)) {
    result <- result & (rowtoMatch[[colName]] == DF[[colName]])
  }

  return(result)

}



#' @title WhichAnswerOriginal
#'
#' @description
#' Given the answers of the students gathered in a table, and a full answer sheet of all versions (Including a "reference/original" version), it finds where those answers are found in the original exam, by copying from the original version the matching rows and binding them in order for every student. It then combines all of them in a list, and includes as well all the remaining student information in the attribute "StudentInfo".
#'
#' It is intended as an internal function to generate the grades, and to identify in a very general way where the answers of the students are (relative to the reference/original version).
#'
#'
#' @param StudentAnswers
#'        DataFrame, each row is a student, each column is some information about said student. Any column not included in \code{names.StudentAnswerQCols} will be understood as information of the student and will be saved as part of the information table when we output the result.
#' @param FullExamAnswerSheet
#'        Answer sheet of all the exam versions, following the conventions of the \code{FullAnswerSheet} outputted by \code{\link{CreateRandomExams}}
#'
#' @param OriginalExamVersion
#'        The version of the original exam, without randomization, as stored on the \code{FullExamAnswerSheet}. The default value is \code{0}, as that is the convention on \code{\link{CreateRandomExams}}
#' @param names.FullExamVersion
#'        The name of the column in which the version of the exam is stored on the \code{FullExamAnswerSheet}. The default value is "\code{Version}", as that is the convention on \code{\link{CreateRandomExams}}
#' @param names.FullExamOriginalCols
#'        The names of the columns that contain the information of the items relative to where they were positioned in the original ordering of the exam, before randomizing the exam. The convention from \code{\link{CreateRandomExams}} is to finish all of them by  "\code{_original}".
#' @param names.CorrectAndIncorrectCols
#'        It should be a character vector. The names of the columns in the \code{FullExamAnswerSheet} that contain the correct and incorrect answers, in that order. This column should have an integer value if it is indeed a correct value in the correct column and a incorrect value in the incorrect value, and \code{NA} otherwise. (The should be "complementary")
#' @param names.StudentAnswerQCols
#'        The names in the \code{StudentAnswers} that store the answers from every student to the exam, ordered. These columns should contain integers values. Where 1 refers to the first answer, and n refers to the nth answers in \strong{their exam}.
#'
#' @param names.StudentAnswerExamVersion
#'        The name of the column in the \code{StudentAnswers} that identifies the version of the exam
#'
#' @details
#' The \code{StudentAnswers} should be a data frame with one student answers represented by every row. The answers of the student to the exam should be ordered.
#'
#'
#' It is important that the colums named \code{names.StudentAnswerQCols} should contain all their answers, if a student didn't answer a question leave a \code{NA} or an invalid integer value as an answer, like 0, or a number larger than the number of answers to that question, so that is is found as out of bounds.
#'
#'
#' @section  Underlying algorithm:
#' To identify the rows on the original exam it does the following:
#'
#' \enumerate{
#' \item It first finds their exam in the full answer sheet by their exam version.
#'
#' \item After that, it  removes from their exam the rows that identify the correct/incorrect choices.
#'
#' \item By trying to match that row with a row on the reference exam it can tell where that quesiton is found on the original exam.
#'
#' \item  Then it identifies where that question is found on the original version, and it finds there which of the possible correct/incorrect choices is found.
#'
#' \item If it didn't find any correct/incorrect choice matching the value given by the student, it marks it as out of bounds and replaces both correct and incorrect columns with \code{NA}.
#'
#' \item If it still doesn't find the row, it simply ignores it, and the output will have one less row.
#'
#' \item Now you can tell how many questions the student answered correctly by looking at how many values are not NA in the correct choice column of the output list.
#'  }
#'
#'
#' @section Removing Questions from the exam:
#'
#' Note that if after creating the exam, you found that a question is bugged and can't be used to grade the exam, all you have to do is tell the student to answer "something" and you only have to remove it from the original/reference version in the Full Answer Sheet. When you apply the grading function, that question will then be ignored.
#'
#' Notice how this creates output lists with different lengths in the case that two students didn't have that same question in their exam.
#'
#' For example, if a exam has 15 questions out of a 50 question document. If student A has a bugged question and student B doesn't, the answer sheet produced for student A will have 14 rows while the one for student B will have 15 rows.
#'
#' @return
#' It returns a list. Each element of the list is a dataframe, and there is one dataframe for each student in the \code{StudentInfo} table provided.
#'
#' All the columns that are not in the columns \code{names.StudentAnswerQCols}  are regarded as "\code{StudentInfo}", and they are added to the attribute "\code{StudentInfo}" of the output as a data frame.
#'
#' \describe{
#'     \item{List elements:}{
#'          They are outputted in order, that is to say, for \code{StudentAnswers[i,]} the list that provides the information for that row will be \code{outputlist[[i]]}.
#'
#'           \code{outputlist[[i]]} is a dataframe that identifies the rows that the student answered as they are found on the original/reference version. Therefore, if a student answeres a certain value, and that value is not reflected on the original version, it get's ignored.
#'
#'     }
#'     \item{\code{StudentInfo} attribute}{A dataframe containing all the student information that wasn't their answers.}
#' }
#'
#'
#' @section Notes:
#' \strong{Note1:} Remember that in the original answer sheet there are two columns, one with correctchoice, another one with wrong choice. If the value is NA of one of those two columns it SHOULD NOT be NA on the other row.
#'
#' \strong{Note2:}  The idea is that the data frames can be read to know the score of the student by counting the number of values that are not NAs on the correct choice column. (The numbers on the correct/incorrect columns themselves can be used for statistical purposes, to tell how many students answered each question).
#'
#' \strong{Note3:} The data frames can be used for many other statistical purposes very easily.
#' @seealso \code{\link{GradeExams}} and \code{\link{ObtainExamStats}} for examples on how to use the output of this function to obtain more detailed information.
#' @export
#' @family Grading exams core Grading Exams
#' @examples
#'
#'
#' asheet_file <-
#'     system.file(
#'         "extdata",
#'         "ExampleTables",
#'         "ExampleAnswerSheet.csv",
#'         package = "TexExamRandomizer")
#' responses_file <-
#'     system.file(
#'         "extdata",
#'         "ExampleTables",
#'         "ExampleResponses.csv",
#'         package = "TexExamRandomizer")
#' FullAnswerSheet <-
#'     read.csv(
#'         asheet_file,
#'         header = TRUE,
#'         stringsAsFactors = FALSE,
#'         na.strings = c("", "NA", "Na"),
#'         strip.white = TRUE)
#' Responses <- read.csv(
#'     responses_file,
#'     header = TRUE,
#'     stringsAsFactors = FALSE,
#'     na.strings = c("", "NA", "Na"),
#'     strip.white = TRUE)
#' compiledanswers <-
#'     WhichAnswerOriginal(
#'         StudentAnswers = Responses,
#'         FullExamAnswerSheet = FullAnswerSheet,
#'         names.StudentAnswerQCols = grep(
#'             names(Responses),
#'             pattern = "^Q.*[[:digit:]]",
#'             value = TRUE),
#'         names.StudentAnswerExamVersion = grep(
#'             names(Responses),
#'             pattern = "Version",
#'             value = TRUE),
#'         OriginalExamVersion = 0,
#'         names.FullExamVersion = "Version",
#'         names.FullExamOriginalCols = grep(
#'             names(FullAnswerSheet),
#'             pattern = "_original",
#'              value = TRUE),
#'         names.CorrectAndIncorrectCols = c(
#'             "choice",
#'             "CorrectChoice")
#'     )
#' nicknames <- attr(compiledanswers, "StudentInfo")$Nickname
#'
#' for (i in 1:length(compiledanswers)) {
#'     cat("Student\t", nicknames[i], " got\t",
#'         sum(!is.na(compiledanswers[[i]]$CorrectChoice)),
#'         " questions correctly\n", sep = "")
#' }

WhichAnswerOriginal <- function(
  StudentAnswers,
  FullExamAnswerSheet,
  OriginalExamVersion = 0,
  names.FullExamVersion = "Version",
  names.FullExamOriginalCols, #Should include last ChoiceCol
  names.CorrectAndIncorrectCols,
  names.StudentAnswerQCols,
  names.StudentAnswerExamVersion) {


  #
  # Small wrapper for the ExtractOriginalRows warning checking:
  #
  # This function only makes sure the int.vector has one element, if it has less than one element it returns 0 and throws a warning,
  # if it has more than one element it returns the first element and throws a warning
  #
  matchSingleElement <- function(int.vector, suppressWarningMessage = FALSE){

    if (length(int.vector) == 0L) {

      if (!suppressWarningMessage) {
        warning("Answer not found in the Answer Sheet")
      }
      return(0L)
    }

    if (length(int.vector) > 1L ) {
      if (!suppressWarningMessage) {
        warning("More than one answer was found matching in the answer sheet\n Returning only the first match")
      }
      return(int.vector[1L])

    }

    return(int.vector)

  }

  # This function should extract in order a data frame from the Original Exam, extracting one row for each of the questions in his answersheet. trying to match what the student answered with the corresponding column in the answer sheet
  # If a student didn't answer a question, the output rows will have that row missing.
  #
  #
  # Args: (Some of the inputs are not given directly since they are given to the outer function and this function is always inside and can understand at all points what those arguments values will be)
  #   OriginalExam: Full Answer sheet for Exam "0", the original exam, before trimming and randomization, the one we want to match our answers against
  #   thisStudent.AnswerSheet: Full answer sheet for the student exam.
  #   names.CorrectAndIncorrectCols: Names of the columns that have the value of that row choice or incorrect choice. (NA if that row is not incorrect in the incorrect choice column and viceversa).
  #
  #
  #   names.FullExamOriginalCols: Character vector, names of Columns that identify the "original" indexes.  They shold allow within one version to identify uniquely one and only one row. The last of these identifiers on the vector should be the identifier that will differentiate between the different choices of that question
  #   thisStudent.Answers: the student answers, should be only one row, all other rows will be ignored, and one column for each question sequentially in order
  #
  # Returns:
  #   A set of rows from OriginalExam ordered in the appearance order of thisStudent.AnswerSheet. Outputting only those rows that are the corresponding rows to the answers in thisStudent.Answers
  #   If the row is not found, it simply doesn't return that row
  ExtractOriginalRows <- function(
    OriginalExam,
    thisStudent.AnswerSheet,
    thisStudent.Answers
  ) {


    thisStudent.UniqueOriginalColumns <- unique(
      thisStudent.AnswerSheet[
        names.FullExamOriginalCols[-length(names.FullExamOriginalCols)]]
      )
    #The last one of the names is suppose to designate the last possible choice, which conserves the original coice
    thisStudent.CombinedLastChoice <- integer(nrow(thisStudent.AnswerSheet))

    # In the mainAnswerSheet there can be more than one Tag that could be matched and divided, NAs should be identified properly
    #
    # CAREFULL TO NOT CHANGE THE NAs OR YOU WILL BREAK THE FUNCTION, SINCE OVER HERE WE ARE JUST ADDING ELEMENTS

    for (colName in names.CorrectAndIncorrectCols) {
      thisStudent.CombinedLastChoice <-
        thisStudent.CombinedLastChoice +
        replace(
          thisStudent.AnswerSheet[[colName]],
          is.na(thisStudent.AnswerSheet[[colName]]),
          0L
        )
    }
    thisStudent.AnswerSheet$CombinedLastChoice <-
      as.integer(
        thisStudent.CombinedLastChoice
      )
    outputIndex <- integer(nrow(thisStudent.UniqueOriginalColumns))

    for (j in seq_along(outputIndex)) {

      #First identify where we can find this answer in our answer sheet
      thisStudent.InfoRowNumber <-
        which(
          FindMatchingRow(
            structure(
              cbind(
                thisStudent.UniqueOriginalColumns[j,],
                thisStudent.Answers[1,j]
              ),
              names = c(names(thisStudent.UniqueOriginalColumns), "CombinedLastChoice")
            ),
            thisStudent.AnswerSheet
          )
        )

      # Then identify that same row in the original exam
      # (or if he doesn't answer at all)
      Original.InfoRowNumber <-
        which(
          FindMatchingRow(
            thisStudent.AnswerSheet[
              matchSingleElement(thisStudent.InfoRowNumber),
              names.FullExamOriginalCols
              ],
            OriginalExam
          )
        )
      outputIndex[j] <- matchSingleElement(Original.InfoRowNumber)
    }

    isOutOfBoundsIndex <- logical(length(outputIndex))
    # If we didn't get a match with the previous method.
    # Maybe it wasn't found because the student wrote an answer which is not taken into account in the answer sheet, but the question does exist.
    # In such cases, we can try to detect it, and mark it as incorrect, by marking as NA both the correct and incorrect columns

    for (j in seq_along(outputIndex)) {
      if (outputIndex[j] == 0) {
        warning("Out Of Bounds answer by the current student")


        #First identify where we can find this answer in our answer sheet
        thisStudent.InfoRowNumber <-
          which(
            FindMatchingRow(
              thisStudent.UniqueOriginalColumns[j,],
              thisStudent.AnswerSheet
            )
          )

        # Then identify that same row in the original exam
        Original.InfoRowNumber <-
          which(
            FindMatchingRow(
              thisStudent.AnswerSheet[
                matchSingleElement(thisStudent.InfoRowNumber, suppressWarningMessage = TRUE),
                names.FullExamOriginalCols
                ],
              OriginalExam
            )
          )
        outputIndex[j] <- matchSingleElement(Original.InfoRowNumber,  suppressWarningMessage = TRUE)
        if (outputIndex[j] != 0L) {
          # If after trying again without the combined last column we actually got a number that it wasn't zero. That means that indeed the answer was out of bounds, but we were able to find the question. Therefore, mark it as out of bounds
          isOutOfBoundsIndex[j] <- TRUE
        }


      }
    }

    TableOfOriginalAnswers <- OriginalExam[outputIndex,]

    # And now we have to check whether there were outofbound indexes to remove the correct and incorrect answers, marking both as NA, since the answer was out of bounds.
    i <- 0
    for (j in seq_along(outputIndex)) {
      if (outputIndex[j] != 0) {
        # Those indices that are not 0 are going to advance one on the table.
        i <- i + 1

        if (isOutOfBoundsIndex[j]) {
          # The convention if the answer was outofbounds is that the answer was actually incorrect, since the student must have written an incorrect answer, that wasn't even in the options given (or he didn't write an answer entirely)
          TableOfOriginalAnswers[i,][names.CorrectAndIncorrectCols] <- NA
        }
      }

    }

    return(TableOfOriginalAnswers)
  }


  StudentAnswersNames <- names(StudentAnswers)
  StudentInfo <- StudentAnswers[, !(StudentAnswersNames %in% names.StudentAnswerQCols)]
  StudentRawAnswers <- StudentAnswers[, names.StudentAnswerQCols]
  StudentAnswerVersionNumber <- StudentAnswers[[names.StudentAnswerExamVersion]]
  OriginalExam <-
    FullExamAnswerSheet[
      FullExamAnswerSheet[[names.FullExamVersion]] == OriginalExamVersion, #rows
      #cols, all of them
      ]

  if (nrow(OriginalExam) == 0) {
    stop("Couldn't find the original exam on the answer sheets")
  }


  forEachRowinStudentAnswers <- function(i){
    #Simply a wrapper to call lapply on the outer function only, rather than concatenating lists
    # Returns:
    #    returns the data frame for one student as described in the documentation of the main function
    thisStudent.ExamVersion <- StudentAnswerVersionNumber[i]
    thisStudent.Answers <- StudentRawAnswers[i,]
    thisStudent.Info <- StudentInfo[i, ]
    thisexam.whichRowMain <-
      FullExamAnswerSheet[[names.FullExamVersion]] == thisStudent.ExamVersion

    if (!(T %in% thisexam.whichRowMain)) {
      # If no match is found:
      warning(
        paste(
          "Version number ",
          thisStudent.ExamVersion,
          " was not found on the original answer sheet.\n",
          "Are you sure you wrote the column names correctly?\n",
          "Student Exam Version column Name: ",
          names.StudentAnswerExamVersion,
          "\n",
          "Original Exam Version column Name: ",
          names.FullExamVersion,
          "\n",
          sep = ""
          )
        )
     # If it didn't find it, return an empty list with zero elements. but the columns of OriginalExam
      return(OriginalExam[0,])
    }

    thisStudent.AnswerSheet <-
      FullExamAnswerSheet[
        thisexam.whichRowMain #row
        , #cols, all of them
        ]


    return(
      ExtractOriginalRows(
        OriginalExam = OriginalExam,
        thisStudent.AnswerSheet = thisStudent.AnswerSheet,
        thisStudent.Answers = thisStudent.Answers
      )
    )
  }

  return(
    structure(
      lapply(1:nrow(StudentAnswers), forEachRowinStudentAnswers),
      StudentInfo = StudentInfo
    )
  )


}

#### GRADING ####

#' @title GradeExams
#'
#' @description
#' Grades an exam given a parsed list by \code{\link{WhichAnswerOriginal}}
#'
#' @param ExamAnswerParsedList List parsed by \code{\link{WhichAnswerOriginal}}
#' @param name.ColCorrect,name.ColIncorrect The names of the correct and incorrect columns in each answer sheet of the \code{ExamAnswerParsedList} respectively.
#' @param MaxOutputGrade Maximum score that one should get if you get a perfect score, before couning the \code{ExtraPoints}
#' @param ExtraPoints Extra points to be added after scoring the exam. This points are added after the scaling is done with \code{MaxOutputGrade}.
#' @param ExtraPointsForAll Scalar numeric value, extra points to be given to all student.
#' @inheritSection WhichAnswerOriginal Removing Questions from the exam
#'
#' @details
#' The score is first added on the base of the number of questions that are found on every parsed list.
#'
#' If a question is removed from an exam, not all students may have that question as explained in the "Removing questions from the exam" section. If the total rows of a certain student list is \eqn{n}, the score is \deqn{c / n * MaxOutputGrade}, where \eqn{c} is the number of correct answers.
#'
#' After that is done, the \code{ExtraPoints} are added.
#' @section Extra Points:
#'
#' The structure of \code{ExtraPoints} and the convention on how the score is calculated taking it into account is worth mentioning in it's own section.

#' The score is calculated as:
#'
#' \deqn{total_{grade} = (c + extra_{all}) / (maxn + extra_{all}) * MaxOutputGrade + extra_{individual}}
#' Where \describe{
#' \item{\code{c}}{Number of correct questions}
#' \item{\code{extra_all}}{Number of extra points for all.
#'
#' This is thought of to be used as a question that you removed from the exam last minute,
#' but that you want to actually count it as correct for every single student. I.e., a question that everyone got correct but it is not taken into consideration in the grading.}
#' \item{\code{extra_individual}}{Number of extra points for that student.}
#' \item{\code{max_n}}{Maximum number of questions in the students exam, which may differ from other students if you had to removed a bugged questions that not everyone had}
#' \item{\code{MaxOutputGrade}}{The scaling to be done. This should be the maximum grade any student "should" get. (The individual extra points are added after the scaling is done)}
#' }
#'
#'
#' @return
#'  It returns the \code{StudentInfo} attribute of the parsed list adding the following columns to it
#'
#'  \describe{
#'      \item{\code{$addedPoints}}{Individual part of ExtraPoints}
#'      \item{\code{$addedAllPoints}}{Extra Points For All}
#'      \item{\code{$maxGrade}}{ Max number of questions for the exam. (It would be different if when removing a question, some students didn't have a question in that exam)}
#'      \item{\code{$Grade}}{Number of correct answers that a student wrote in an exam}
#'      \item{\code{$Grade_Total_Exam}}{This is the \code{total_grade} as explained on the Extra Points section.}
#'  }
#' @export
#' @family Grading Exams
#' @examples
#'
#' #First part coming from FindMatchingRow example
#'
#' asheet_file <-
#'     system.file(
#'         "extdata",
#'         "ExampleTables",
#'         "ExampleAnswerSheet.csv",
#'         package = "TexExamRandomizer")
#' responses_file <-
#'     system.file(
#'         "extdata",
#'         "ExampleTables",
#'         "ExampleResponses.csv",
#'         package = "TexExamRandomizer")
#' FullAnswerSheet <-
#'     read.csv(
#'         asheet_file,
#'         header = TRUE,
#'         stringsAsFactors = FALSE,
#'         na.strings = c("", "NA", "Na"),
#'         strip.white = TRUE)
#' Responses <- read.csv(
#'     responses_file,
#'     header = TRUE,
#'     stringsAsFactors = FALSE,
#'     na.strings = c("", "NA", "Na"),
#'     strip.white = TRUE)
#' compiledanswers <-
#'     WhichAnswerOriginal(
#'         StudentAnswers = Responses,
#'         FullExamAnswerSheet = FullAnswerSheet,
#'         names.StudentAnswerQCols = grep(
#'             names(Responses),
#'             pattern = "^Q.*[[:digit:]]",
#'             value = TRUE),
#'         names.StudentAnswerExamVersion = grep(
#'             names(Responses),
#'             pattern = "Version",
#'             value = TRUE),
#'         OriginalExamVersion = 0,
#'         names.FullExamVersion = "Version",
#'         names.FullExamOriginalCols = grep(
#'             names(FullAnswerSheet),
#'             pattern = "_original",
#'             value = TRUE),
#'         names.CorrectAndIncorrectCols = c(
#'             "choice",
#'             "CorrectChoice")
#'     )
#' # Actual Code
#'
#'
#' ExtraPoints_individual <- runif(nrow(Responses), min = 1, max = 10)
#' ExtraPoints_forall <- 2
#' GradedStudentTable <-
#'     GradeExams(
#'         compiledanswers,
#'         name.ColCorrect = "CorrectChoice",
#'         name.ColIncorrect = "choice",
#'         MaxOutputGrade = 100,
#'         ExtraPoints = ExtraPoints_individual,
#'         ExtraPointsForAll = ExtraPoints_forall
#'     )
#'
#'
#'
#'
GradeExams <- function(
  ExamAnswerParsedList,
  name.ColCorrect,
  name.ColIncorrect,
  MaxOutputGrade = 100,
  ExtraPoints = 0,
  ExtraPointsForAll = 0) {


  GradeDF <- function(DataFrame) {

    # Thanks to the formatting of the original answer sheet, we only need to count when there is not an NA in the original data frame.

    # TODO: Improve grading to add question values? This would need to be added to the core of the exam.
    return(
      structure(
        sum(
          !is.na(
            DataFrame[[name.ColCorrect]]
          )
        ),
        MaxGrade = nrow(DataFrame)
      )
    )
  }

  assertthat::assert_that(
    length(ExamAnswerParsedList) == nrow(attr(ExamAnswerParsedList, "StudentInfo"))
  ) #If not, this list wasn't parsed from one of our functions...

  Grade <- numeric(length(ExamAnswerParsedList))
  maxGrade <- numeric(length(ExamAnswerParsedList))
  for (i in seq_along(ExamAnswerParsedList)) {
    tempGrade <-  GradeDF(ExamAnswerParsedList[[i]])
    Grade[i] <- tempGrade
    maxGrade[i] <- attr(tempGrade, "MaxGrade")

  }

  output <- attr(ExamAnswerParsedList, "StudentInfo")



  if (is.null(ExtraPoints)) {
      IndividualExtraPoints <- rep(0, length(ExamAnswerParsedList))
  } else if  (length(ExtraPoints) == 1) {
      IndividualExtraPoints <- rep(ExtraPoints, length(ExamAnswerParsedList))
  } else if (length(ExtraPoints) != length(ExamAnswerParsedList)) {
      stop("ExtraPoints should have the same length as the number of students passed in.")
  } else {
      IndividualExtraPoints <- ExtraPoints
  }

  if (length(ExtraPointsForAll) > 1) {
      warning("ExtraPointsForAll should be a scalar, selecting only the first element")
      ExtraPointsForAll <- ExtraPointsForAll[1]
  } else if (is.null(ExtraPointsForAll)) {
      ExtraPointsForAll <- 0L
  }

  ExtraPointsForAll <- rep(ExtraPointsForAll[1], length(ExamAnswerParsedList))


  if (any(0 %in% maxGrade)) {
    warning("Some answer sheets have no answers")
    maxGrade <- replace(x = maxGrade, list = 0, values = 1)
  }
  output$addedPoints <- IndividualExtraPoints
  output$addedAllPoints <- ExtraPointsForAll
  output$maxGrade <- maxGrade

  output$Grade <- Grade
  output$Grade_Total_Exam <-
    (output$Grade + output$addedAllPoints) / (output$maxGrade + output$addedAllPoints) * MaxOutputGrade + output$addedPoints

  return(output)


}

#### STATS ####


#' Obtaining exam statistics
#'
#' This function gets an answer sheet of the original version of the exam as a data frame, and a parsed list, which is obtained from \code{\link{GradeExams}} and it outputs the statistics of how many answers are parsed exam, that is graded and obtains from there
#'
#' @param OriginalExamAnswerSheet The answer sheet of the original exam. (In this package the convention is the exam version "0")
#' @param ExamAnswerParsedList A parsed list for every student, as outputted by \code{\link{GradeExams}}
#' @param names.FullExamOriginalCols Names of those columns that in the answer sheet identify for all versions where that item is found on the original columns, (i.e., as ordered from the original version exam)
#'
#'
#' @return Returns the \code{OriginalExamAnswerSheet} with a column added to it, named "\code{ExamAnswerCount}" that counts the number of answers for each question
#'
#'
#' @export
#' @family Grading Exams
#' @examples
#'
#' asheet_file <-
#'     system.file(
#'         "extdata",
#'         "ExampleTables",
#'         "ExampleAnswerSheet.csv",
#'         package = "TexExamRandomizer")
#' responses_file <-
#'     system.file(
#'         "extdata",
#'         "ExampleTables",
#'         "ExampleResponses.csv",
#'         package = "TexExamRandomizer")
#' FullAnswerSheet <-
#'     read.csv(
#'         asheet_file,
#'         header = TRUE,
#'         stringsAsFactors = FALSE,
#'         na.strings = c("", "NA", "Na"),
#'         strip.white = TRUE)
#' Responses <- read.csv(
#'     responses_file,
#'     header = TRUE,
#'     stringsAsFactors = FALSE,
#'     na.strings = c("", "NA", "Na"),
#'     strip.white = TRUE)
#' compiledanswers <-
#'     WhichAnswerOriginal(
#'         StudentAnswers = Responses,
#'         FullExamAnswerSheet = FullAnswerSheet,
#'         names.StudentAnswerQCols = grep(
#'             names(Responses),
#'             pattern = "^Q.*[[:digit:]]",
#'             value = TRUE),
#'         names.StudentAnswerExamVersion = grep(
#'             names(Responses),
#'             pattern = "Version",
#'             value = TRUE),
#'         OriginalExamVersion = 0,
#'         names.FullExamVersion = "Version",
#'         names.FullExamOriginalCols = grep(
#'             names(FullAnswerSheet),
#'             pattern = "_original",
#'             value = TRUE),
#'         names.CorrectAndIncorrectCols = c(
#'             "choice",
#'             "CorrectChoice")
#'     )
#' OriginalAnswerSheet <- FullAnswerSheet[FullAnswerSheet$Version == 0,]
#' ExamStats <-
#'     ObtainExamStats(
#'         OriginalExamAnswerSheet = OriginalAnswerSheet,
#'         ExamAnswerParsedList = compiledanswers,
#'         names.FullExamOriginalCols =  grep(
#'             names(FullAnswerSheet),
#'             pattern = "_original",
#'             value = TRUE)
#'     )


ObtainExamStats <- function(
  OriginalExamAnswerSheet,
  ExamAnswerParsedList,
  names.FullExamOriginalCols
){

  ExamAnswerCount <- integer(nrow(OriginalExamAnswerSheet))
  # Defining it before AddQuestionCountFromDF to prevent warnings form syntax Rstudio thingy

  AddQuestionCountFromDF <- function(AnswerSheet){
    # this function takes the answer sheet and updates the ExamAnswerCount of the main function.
    # Args:
    #   AnswerSheet, the answer sheet form which to update the count
    # Returns:
    #   NULL if the answer sheet is empty, 1L if it succeded in complete the function
    if (nrow(AnswerSheet) == 0) {
      return(NULL)
    }
    for (i in 1:nrow(AnswerSheet)) {
      tempfind.logical <- FindMatchingRow(AnswerSheet[i,names.FullExamOriginalCols], OriginalExamAnswerSheet)
      if (length(tempfind.logical) == length(ExamAnswerCount)) {
        ExamAnswerCount <<-
          ExamAnswerCount +
          as.integer(
            FindMatchingRow(AnswerSheet[i,names.FullExamOriginalCols], OriginalExamAnswerSheet)
          )
      } else if (length(tempfind.logical) == 0) {
        warning("FindMatchingRow didn't match any row")
      }
    }

    return(1L)

  }


  for (parsedStudentAnswer in ExamAnswerParsedList) {

    AddQuestionCountFromDF(parsedStudentAnswer)
    #This function updates internally ExamAnswerCount to count the new answersheet form a student
  }

  return(
    cbind(
      OriginalExamAnswerSheet,
      ExamAnswerCount = ExamAnswerCount
    )
  )


}
