# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017


#' @title ReplacePreambleCommand
#' @description
#' This functions gets a character vector in which each element represents a line of
#' a preamble of a 'LaTeX' document, and it replaces the definition of the command \code{\\commandName} to have the value \code{commandValue}.
#' @details It only modifies the value of the command by replacing instances of
#'
#' \code{\\newcommand\{\\commandName\}\{<previous definition>\}} with instances of
#'
#' \code{\\newcommand\{\\commandName\}\{<commandValue>\}}.
#'
#' @details
#' Keep in mind that both \code{commandName} and \code{commandValue} are placed directly inside a regex.
#'
#' If you want to "hide" a certain definition of a command from being found and replaced by this function,  simply define it by using \code{\\def} or \code{\\newcommand*} or a \code{\\renewcommand} when you define them.
#'
#'
#' Make sure you are using a one-line definition in commands that you want replaced, since this won't be able to detect commands that are defined in multiple lines in 'LaTeX'.
#'
#' Also, note how certain invalid things in 'LaTeX' would still be matched by this regex,
#' however you should find those errors before you start using this program since those errors would not allow you to compile the 'LaTeX' document on the first place.
#'
#' Lastly, if it doesn't find a command on the document, it silently ignores it.
#'
#' @param x
#'       A character vector, each element is suppose to represent a line
#' @param commandName
#'       A string identifying either the command name
#' @param commandValue
#'       Replacement for the definition of commandName
#' @return A character vector, with the preamble, replacing all instances of
#' \code{ \\newcommand\\commandName\{<random text>\}} with
#' \code{ \\newcommand\\commandName\{commandValue\}}
#' @export
#' @family Preamble adjustment
#' @examples
#' new_preamble <- ReplacePreambleCommand( TexExamRandomizer::testdoc$preamble, "nickname", "Alex")

ReplacePreambleCommand <- function(x, commandName, commandValue){



  # Intention when using the balanced thing.
  # (\\{(?:[^{}%]|(?1))*\\})\z
  #


  pattern1 <- "^([^%]*)\\\\newcommand(\\{)?\\\\"
  pattern2 <- "(\\})?(\\{(?:[^{}%]|(?4))*\\})" #The second part matches balanced parenthesis


  pattern <- paste(pattern1, commandName, pattern2, sep = "")

  replacement <- paste("\\1\\\\newcommand{\\\\", commandName, "}{", commandValue, "}", sep = "")

    return(
      sub(pattern = pattern,
          replacement = replacement,
          x = x,
          perl = TRUE)
    )
}


#' @title ReplaceFromTable
#' @description
#'
#' Given a 'LaTeX' file represented as a character vecotr with \code{x}, it replaces from a table the commands given by \code{commandNames}. for the values found on the table.
#'
#'
#' \code{\\newcommand\{\\commandName[i]\}\{table[tableRow, columnName[i]]\}}.
#'
#' @details To do the replacement for each item, it uses the function \code{\link{ReplacePreambleCommand}}. See the details in that function for more information.
#'
#' @param x
#'       A character vector, each element is suppose to represent a line
#' @param table Data frame from which to extract the information
#' @param tableRow Integer, row of the \code{table} to be used
#' @param columnNames Character vector with the names of the columns to be used
#' @param commandNames Character vector with the same length as \code{columnNames}. Contains the names of the 'LaTeX' commands to be replaced.
#' @return A character vector, representing the text \code{x}, where all instances of
#' \code{ \\newcommand\\commandNames[i]\{<random text>\}} have been replaced with
#' \code{ \\newcommand\\commandNames[i]\{table[tableRow, columnName[i]\}}.
#'
#' @export
#' @family Preamble adjustment
#' @examples
#'
#' custom_preambles <- list()
#' for (i in 1:nrow(TexExamRandomizer::testclass)) {
#'     custom_preambles <-
#'         c(
#'             custom_preambles,
#'             list(
#'                 TexExamRandomizer::ReplaceFromTable(
#'                     TexExamRandomizer::testdoc$preamble,
#'                     table = TexExamRandomizer::testclass,
#'                     tableRow = i,
#'                     columnNames = c("Class", "Roll.Number", "Nickname"),
#'                     commandNames = c("class", "rollnumber", "nickname")
#'                 )
#'             )
#'
#'         )
#'
#' }

ReplaceFromTable <- function(x, table, tableRow, columnNames, commandNames){

  # This function  applies ReplacePreambleCommand to a series of commands, which information we find  in a row of a data frame.
  # Args:
  #   x: Character vector in which to replace things
  #   table: data frame from which to extract the information
  #   tableRow: Row to extract from
  #   commandNames: Command names defined as  \newcommand{\commandName}{ bla bla bfadsfads} to be replaced
  #   columnNames: Column from which to replace the value of the command. It will end up looking as
  #     \newcommand{\commandName}{<table[tableRow, columnName[i]]>}
  #
  # Returns:
  #   Returns x with the values of the commands replaced by their new values found on the table


  assertthat::assert_that(length(commandNames) == length(columnNames))

  for (i in seq_along(columnNames)) {
    # x %<>% ReplacePreambleCommand(commandNames[i], table[[columnNames[i]]][tableRow])
    x <-
      ReplacePreambleCommand(
        x,
        commandNames[i],
        table[[columnNames[i]]][tableRow]
      )
  }

  return(x)
}
