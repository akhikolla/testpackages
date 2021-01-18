# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017

# Functions required to create the Structure

#' @title IsWellSectioned
#' @description Function to assure a set of sections is well sectioned.
#'
#' @details Basically it makes sure that, \eqn{u[1]<v[1]<u[2]<v[2]}, etc
#'
#' @param u Vector, it assumes it is ordered in ascending ordered
#' @param v Vector, it assumes it is ordered in ascending ordered
#
#' @return Logical value, TRUE if it is well ordered, FALSE it is not
#'
#' @author Alejandro Recuenco \email{alejandrogonzalezrecuenco@@gmail.com}
#' @keywords internal
#' @family Structuring Document
IsWellSectioned <- function(u, v) {


    i <- 1
    n <- length(u)
    if (n != length(v)) {
        return(FALSE)
    }

    if (n == 0) {
        warning("The sectioning lines are empty")
        return(TRUE)
    }
    first_num  <- u[1]
    second_num <- v[1]
    nexttocheckisv <- FALSE
    while (i <= n) {

        if (second_num < first_num) {
            return(FALSE)

        } else {
            first_num <- second_num
            if (!nexttocheckisv) {
                if (i == n) {
                    return(TRUE)
                }
                i <- i + 1
                second_num <- u[i]

            } else {

                second_num <- v[i]
            }
            nexttocheckisv <- !nexttocheckisv

        }

    }

    return(TRUE)

}



#' @title Structuring functions
#' @name FindStructure
#' @description These internal functions provide functionality to find environment names and their use
#' @param x string vector, each line should represent one line of a text file or a section fo a text file.
#' @param cmdName Command to search for in \code{x}.
#'
#' @details In \code{\link{FindBegin}} and \code{\link{FindEnd}}, \code{cmdName} refers to the name of the command that would start an environment. Following the 'LaTeX' convention of "\\begin\{\code{cmdName}\}" or "\\end\{\code{cmdName}\}" respectively. However, it is not a full throrough check.
#' They will only be found by this class if they are only preceded by alphanumeric characters and spaces, this is to force the user to use begin and end environments at the start of a new line. TODO: CONSIDER AYBE CHANGING THIS THIS, IT IS EASY.
#'
#' On the other hand, in the function  \code{\link{FindCommand}}. It finds the command "\\\code{cmdName}". And in this case it is less rescrittive, as long as the line is not commented, it will find it.
#' Make sure to not write slashes before the "\\\code{cmdName}", since you might bug the program if it thinks you wrote the command but you just wrote some slashes and then the command name afterwards.
#'
#' All functions don't search for commands if the commands have been commented with the latex comment command in the same line... don't try to use multiple line comments on latex please...
#'
#' \strong{IMPORTANT}, instead of just writing something alphanumeric, these function actually use the  \code{cmdName} as a regular expression, which might be useful in many cases, but be careful with this.
#'
#' TODO: Implement options for regular expressions

#' @family Structuring Document
NULL


#' @inheritParams FindStructure
#' @return Returns a numeric vector, indicating each occurrance of a start of a environment that looks like \\begin\{\code{cmdName}\}
#'
#' @keywords internal
#' @rdname FindStructure
FindBegin <- function(x, cmdName) {

    # Double % because of sprintf removing one
    pattern = "^[^%%]*\\\\begin\\{%s\\}"

    return(
        which(
            grepl(
                pattern = sprintf(pattern, cmdName),
                x = x
            )
        )
    )
}

#'The function\code{\link{FindEnd}} is function returns the position  of the vector string x in which it finds a match to the ending of an latex environment with name \code{cmdName}
#' @inheritParams FindStructure
#' @return Returns a numeric vector, indicating each occurrence of a start of a environment that looks like \\end\{\code{cmdName}\}
#'
#' @keywords internal
#' @rdname FindStructure

FindEnd <- function(x, cmdName) {

    # Double % because of sprintf removing one
    pattern = "^[^%%]*\\\\end\\{%s\\}"

    return(
        which(
            grepl(
                pattern = sprintf(pattern, cmdName),
                x = x
            )
        )
    )

}

#' @inheritParams FindStructure
#' @return Returns a numeric vector, indicating each occurrence of the command \\\code{cmdName} found in the document.
#' @keywords internal
#'
#' @rdname FindStructure

FindCommand <- function(x, cmdName) {

    pattern1 <- "^[^%]*\\\\"
    pattern2 <- "([^[:alpha:]]|$)"

    return(
        which(
            grepl(
                pattern = paste(pattern1, cmdName, pattern2, sep = ""),
                x = x
            )
        )
    )
}


#' @title DivideFile
#' @description Function that takes a vector of text lines, \code{x}, and divides it in preamble and document.
#' @details It ignores everything after the first \\end\{document\} and it will throw and error if it finds more than one \\begin\{document\} before that
#' @param x A character vector, each element represents one line of the latex document
#' @return  Returns a list with two character vectors:
#' \describe{
#' \item{preamble}{ A character vector that includes \emph{every line} of \code{x} up to \\begin\{document\}}
#' \item{document}{A character vector that includes \emph{every line} from \\begin{document} to the first \\end\{document\}} }
#' @export
#' @family Structuring Document
#' @examples
#' file <- system.file(
#'     "extdata",
#'     "ExampleTexDocuments",
#'     "exam_testing_jsonparser.tex",
#'     package = "TexExamRandomizer"
#' )
#' x <- readLines(file)
#' DivideFile(x)
DivideFile <- function(x) {

    beginDocLines <- FindBegin(x, "document")
    endDocLines   <- FindEnd(x, "document")[1]

    beginDocLines <- beginDocLines[beginDocLines < endDocLines]

    assertthat::assert_that(IsWellSectioned(beginDocLines, endDocLines))

    #beginDocLines should not be line 1
    assertthat::assert_that(beginDocLines > 1)

    # preamble <- x[1:(beginDocLines - 1)]
    #
    # document <- x[beginDocLines:endDocLines]
    #
    return(list(preamble = x[1:(beginDocLines - 1)],
                document =  x[beginDocLines:endDocLines]))

}

#' @title Structure Document
#'
#' @description Function that takes a character vector, \code{x}, representing a 'LaTeX' file and it outputs a tree structure with the structure specified by \code{layersNames} and \code{layersCmd}.
#'
#' It assumes \code{x} is representing a 'LaTeX' file that can has been checked it compiles apropitaly before we make anymodification.
#'
#' Note however that this function only moves lines around, it doesn't split a line in two.
#' @inheritParams DivideFile
#' @param layersNames A character vector, with each element representating the environment name to be searched as \code{cmdName} as describe in \code{\link{FindBegin}} and \code{\link{FindEnd}}
#' @param layersCmd A character vector, with the same length as \code{layersNames}. with each element representing the environment command to be serached as \code{cmdName} as described in \code{\link{FindCommand}}.
#'
#'
#' @details
#' Both \code{layersNames} and \code{layersCmd} must have the same length, since for each index, \code{i}, \code{layersNames[i]} and \code{layersCmd[i]} refer to one layer of the tree structure of the document. Consequent layers must be found inside previous layers.
#'
#' If it finds the structure of the document to not be completed, it will throw an error.
#'
#' @return  It returns a list, with each element having a name. Recreating the tree structure identified by  \code{layersNames} and \code{layersCmd} in the text file \code{x}.
#'
#' It first divides the document into two lists:
#' \describe{
#' \item{preamble}{Contains a character vector identifying everything before the \\begin\{document\}}
#' \item{document}{Contains the tree structure identifying the document}
#' }
#'
#' Now, the naming convention for each layer of the document is as follows. We will use the convention \code{<layerName>},  \code{<layerCmd>}.
#'
#'  Note the convention first, everything that it finds prior to the first environment, it throws it into a character vector that it calls \code{prior_to_<layesName>}.
#'  After the first environment \code{<layerName>} ends, it assumes that everything from that \code{\\end\{<layerName>\}} onwards corresponding to the next environment, and it will throw it to the prior part of that one.
#'  \code{post_to_<layerName>}
#' \describe{
#' \item{\code{prior_to_layersName}}{Includes everything up to the first \code{\\begin\{<layerName>} without including that line}
#' \item{\code{1_<layerName>_begin_<layerName>}}{
#'   Includes the \code{\\begin\{layerName\}} for the 1st section, and everything until it finds the first \code{\\<layerCmd>}}
#' \item{\code{1_<layerName>_1_<layerCmd>}}{
#'   Includes everything from the 1\eqn{^{st}} \code{\\<layerCmd>} until the second \code{\\<layerCmd>}, without including the line in which the second command is found
#' }
#' \item{\code{1_<layerName>_2_<layerCmd>}}{
#'  Same thing... and it keeps going until the last  \code{\\<layerCmd>} is found
#' }
#' \item{\code{1_<layerName>_end_<layerName>}}{
#' It includes the \code{\\end\{<layerName>\}} for the 1st section.
#' }
#' \item{...}{
#' It then repeats the same structure for the next environment, changing the naming convention to start with 2_<...> and so on until it does the last environemt}
#' \item{\code{post_to_<layerName>}}{
#' After the last layer ends with \code{\\end\{layerName\}}, it throws the rest of the lines into this last character vector}
#'
#' }
#'
#' This structure is applied recursively to each \code{i_<layerName>_j_<layerCmd>} of the previous layer to find the structure for the next layer.
#' The result is a tree of lists, with names that identify the whole structure, and the ending node of each branch is always a character vector
#'
#' \strong{IMPORTANT NOTE:} Note that this function only rearranges the lines of the document, it can't split a document between a line. So if you want to make sure something always stays together, put them both in the same line. This is intentional, to force a more clear structure on the document that will be parsed
#'
#' In Summary, the sketch of the tree structure would be:
#' \itemize{
#'   \item preamble
#'   \item Document
#'         \itemize{
#'         \item prior_to_LayerName[1]
#'         \item 1_layerName[1]_begin_layerName[1]
#'         \item 1_layerName[1]_1_layerCmd[1]
#'              \itemize{
#'              \item prior_to_LayerName[2]
#'              \item 1_layerName[2]_begin_layerName[2]
#'              \item 1_layerName[2]_1_layerCmd[2]
#'                   \itemize{
#'                   \item Continues...
#'                   }
#'              \item 1_layerName[2]_2_layerCmd[2]
#'                   \itemize{
#'                   \item Continues...
#'                   }
#'              \item ...
#'              \item  post_to_layerName[2]
#'              }
#'         \item 2_layerName[1]_begin_layerName[1]
#'         \item 2_layerName[1]_1_layerCmd[1]
#'              \itemize{
#'              \item ...
#'              }
#'        \item ...
#'        \item n_layerName[1]_end_layerName[1]
#'        \item  post_to_layerName[1]
#'         }
#'
#' }
#'
#' If a  \code{\\<layerCmd>} is not found inside an environment, everything inside that environment is thrown into the begin_layerName part and instead of the numbered environments, an empty character list is added in the middle, with name \code{empty_<layerCmd>} section.
#'
#' @seealso \link{FindStructure} for more information on the details of how the layers are found.
#' @export
#' @family Structuring Document
#' @examples
#' file <- system.file(
#'     "extdata",
#'     "ExampleTexDocuments",
#'     "exam_testing_jsonparser.tex",
#'     package = "TexExamRandomizer"
#' )
#' x <- readLines(file)
#' layersNames <- c("questions", "choices")
#' layersCmd <- c("question", "(choice|CorrectChoice)")
#' doc <- StructureDocument(x, layersNames, layersCmd)
StructureDocument <- function(x, layersNames, layersCmd){

    # Checking
    assertthat::assert_that(!any(grepl("[[:digit:]]", layersNames))) # Making sure the layerNames have no digits
    assertthat::assert_that(!any(grepl("[[:digit:]]", layersCmd))) # Making sure the layersCmd have no digits
    assertthat::assert_that(assertthat::are_equal(length(layersNames), length(layersCmd))) #Same Length

    assertthat::assert_that(
        length(layersNames) == length(unique(layersNames)),
        length(layersCmd) == length(unique(layersCmd)),
        msg = "Aborting, layers must have unique names, you can't have nested layers with the same names"
        )


    # Code of StructureExam

    Document <- DivideFile(x)

    Document[[2]]  <- CompileDocument(Document[[2]], layersNames, layersCmd)


    return(Document)


}
