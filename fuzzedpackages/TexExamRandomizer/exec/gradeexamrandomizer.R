#!/usr/local/bin/Rscript
#### SCRIPT ####
#### Grading exams,
# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017
#### Functions ####



StopIfTableNotFound <- function(filename) {
  # Wrap up function to cause the error calling if the filename is not just a single table
  # Returns invisible T if the function did not throw an error
  if (length(filename) != 1) {
    stop("The filename found is not just one table, there is a conflict")
  }
  if (!grepl(pattern = "\\.csv$", x = filename)) {
    stop("The filename given is not a csv file, consider changing the saving extension")
  }

  invisible(T)

}


#### Definitions assumed during this script

ASHEET_FILE_PATTERN <- "answer[^/\\]*sheet[^/\\]*\\.csv$"
# First the answer sheet is detected, then the responses between all the files remaining in the folder
RESPONSES_FILE_PATTERN <- "(response|answer)[^/\\]*\\.csv$"

ASHEET_COLCORRECT <-  "CorrectChoice"
ASHEET_COLINCORRECT <- "choice"
# Assumes only two columns, although main program can do more than 2 columns
ASHEET_VERSION <- "Version"
ASHEET_VERSION_NUMBER <- 0
ASHEET_ORIGINAL_PATTERN <- "_original$"
RESP_QUESTIONCOL_PATTERN <- "Q.*[[:digit:]]+[^[:alpha:]]*$"
RESP_VERSION_PATTERN <- "version"
RESP_EXTRAPOINTS_PATTERN <- "(?i)Extra.*Points"

MAX_OUTPUT_GRADE <- 100
# For iconthailand
RESP_STUDENTCLASS_PATTERN <- "class"

# require(TexExamRandomizer)
library(assertthat)
library(TexExamRandomizer)
library(optparse)
option_list <- list(
  make_option(
    c("-r", "--resp"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Filename of the responses of the students"
  ),
  make_option(
    c("-a", "--answer"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Filename of the answersheet"
  ),
  make_option(
      c("-d", "--dir"),
      action = "store",
      default = NA,
      type = 'character',
      help = "Filename of a directory in which to find the answers and filename with an specific naming convention"
  ),
  make_option(
      c("--max"),
      action = "store",
      default = NULL,
      type = 'integer',
      help = "Maximum grade scale. The default is 100"
  )
)

opt <-
  parse_args(
    OptionParser(option_list = option_list),
    positional_arguments = T
  )

# Checking main tex file exists

if (is.na(opt$options$dir)) {
  # Not passing directory as an option
  if (is.na(opt$options$resp) || is.na(opt$options$answer)) {
    stop("Exam sheet or answer sheet not specified, neither was the directory")
  } else {
    fileM <- opt$options$resp # Students' responses
    answerSheet <- opt$options$answer # Ansewr sheet
    folder <- dirname(fileM) # Output directory is set to be the place where the respones are found by convention
  }
} else {
  folder <- opt$options$dir
  if (!is.dir(folder)) {
    # IF a file is introduced as input of a folder it can try to use the folder in which the file is found instead...
    folder <- dirname(folder)
  }
  filenames <- dir(folder, full.names = T)
  which.answerSheet <-
    grepl(
      pattern = ASHEET_FILE_PATTERN,
      x = filenames,
      ignore.case = T
    )
  answerSheet <- filenames[which.answerSheet]
  # Remove the answer sheet from the filenames, we wouldn't want to pick it twice, look how the grepl for fileM could also pick the answer sheet if we were not careful
  filenames <- filenames[!which.answerSheet]
  fileM <-
    filenames[
      grepl(
        pattern = RESPONSES_FILE_PATTERN,
        x = filenames,
        ignore.case = T
      )
      ]

}



if (!is.null(opt$options$max)) {
    MAX_OUTPUT_GRADE <- as.integer(opt$options$max)
}


StopIfTableNotFound(fileM)
StopIfTableNotFound(answerSheet)

# Now error checking to make sure the things over ehre are as we intended them
f.StudentSheet <- fileM
f.AnswerSheet <- answerSheet







cat("Reading tables...\n", file = stderr())



StudentSheet <-
  read.table(
    f.StudentSheet,
    header = TRUE, sep = ",", stringsAsFactors = FALSE,
    na.strings = c("","NaN","NA","Na","NAN"), strip.white = TRUE
  )
AnswerSheet <-
  read.table(
    f.AnswerSheet,
    header = TRUE, sep = ",", stringsAsFactors = FALSE,
    na.strings = c("","NaN","NA","Na","NAN"), strip.white = TRUE
  )


colsStudentSheet <- names(StudentSheet)


colsAnswerSheet <- names(AnswerSheet)

names.StudentAnswerQCols <-
  colsStudentSheet[
    grepl(
      pattern = RESP_QUESTIONCOL_PATTERN,
      x = colsStudentSheet,
      ignore.case = T
    )
  ]


names.StudentAnswerExamVersion <-
  colsStudentSheet[
    grepl(
      pattern = RESP_VERSION_PATTERN,
      x = colsStudentSheet,
      ignore.case = T
    )
    ]


ExtraPointsCol <- grep(names(StudentSheet), pattern = RESP_EXTRAPOINTS_PATTERN, perl = TRUE)

if (length(ExtraPointsCol) > 0L) {
    ExtraPoints <- StudentSheet[[ExtraPoints[1]]]
} else {
    ExtraPoints <- 0
}


# TODO: PUT this three checks into just one function wrapper with some text to announce what it is

if (length(names.StudentAnswerQCols) == 0L) {
  stop("Question not found in responses sheet, make sure to write a column with the necessary name, using \"version\"")
}


if (length(names.StudentAnswerExamVersion) == 0L) {
  stop("Exam  version not found in responses sheet, make sure to write a column with the necessary name, using \"version\"")
} else if (length(names.StudentAnswerExamVersion) != 1L) {
  warning( "More than one exam matches the condition for the column, picking the first one")
  names.StudentAnswerExamVersion <- names.StudentAnswerExamVersion[1L]
}



cat("\tCompiling answers...\n", file = stderr())

ListCompiledAnswers <-
  WhichAnswerOriginal(
    StudentAnswers = StudentSheet,
    FullExamAnswerSheet = AnswerSheet,
    OriginalExamVersion = ASHEET_VERSION_NUMBER,
    names.FullExamVersion = ASHEET_VERSION,
    names.FullExamOriginalCols =
      colsAnswerSheet[
        grepl(colsAnswerSheet, pattern = ASHEET_ORIGINAL_PATTERN)
        ],
    names.CorrectAndIncorrectCols =
      tail(
        x = colsAnswerSheet,
        n = 2L
      ),
    names.StudentAnswerQCols = names.StudentAnswerQCols,
    names.StudentAnswerExamVersion = names.StudentAnswerExamVersion
  )


cat("\tGrading exams...\n", file = stderr())

AllAnswers <-
  GradeExams(
    ExamAnswerParsedList = ListCompiledAnswers,
    name.ColCorrect = ASHEET_COLCORRECT,
    name.ColIncorrect = ASHEET_COLINCORRECT,
    ExtraPoints = ExtraPoints,
    MaxOutputGrade = MAX_OUTPUT_GRADE
  )


cat("Generating stats...\n", file = stderr())

ExamStats <-
  ObtainExamStats(
    OriginalExamAnswerSheet =
      AnswerSheet[AnswerSheet$Version == ASHEET_VERSION_NUMBER, ],
    ExamAnswerParsedList = ListCompiledAnswers,
    names.FullExamOriginalCols =
      colsAnswerSheet[
        grepl(colsAnswerSheet, pattern = ASHEET_ORIGINAL_PATTERN)
        ]
  )


baseName <-
  sub(
    pattern = paste(
        "[[:space:]]*",
        RESPONSES_FILE_PATTERN,
        sep = ""
    ),
    replacement = "",
    ignore.case = T,
    x = basename(fileM)
  )

file.AllAnswers <-
  file.path(
    folder,
    paste(baseName, "_Graded.csv", sep = "")
  )

file.ExamStats <-
  file.path(
    folder,
    paste(baseName, "_Stats.csv", sep = "")
  )

cat("Saving grades... \n", file = stderr())

write.table(
  x = AllAnswers,
  file = file.AllAnswers,
  sep = ",",
  quote = FALSE,
  na = "",
  row.names = FALSE
)

cat("Saving stats...\t", file = stderr())

write.table(
    x = ExamStats,
    file = file.ExamStats,
    sep = ",",
    quote = FALSE,
    na = "",
    row.names = FALSE
)


cat("done\n", file = stderr())
