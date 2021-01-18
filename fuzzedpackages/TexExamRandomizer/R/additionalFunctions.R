# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017
# Subsetting
#
#
#
#' @title Subsetting Document
#'
#' @description If you use a function such as A[vector], it only preserves the names attribute, it doesn't preserve any other attribute. This Function attempts to "fix" that manually, for lack of a better solution, instead of using S3 classes. Therefore, it can be used with other documents.
#'
#' It subsets the attributes that have the same length as the Document itself and preserve all other attributes with length 0
#'
#' @param Document Document that we want to subset
#' @param vector The vector that provides the subsetting, like \code{`[`}
#'
#' @keywords internal
#' @details Although this function was created with documents in mind, it can work with any R array that holds attributes and has a subsetting \code{`[`} function call.
#'
SubsetWithAtributes <- function(Document, vector) {

    n <- length(Document)

    result <- Document[vector]

    if (length(result) == 0) {
        return(result)
    }

    attribute_list <- attributes(Document)
    name_attribute <- names(attribute_list)

    for (i in seq_along(attribute_list)) {
        if (length(attribute_list[[i]]) == n) {
            attr(result, name_attribute[i]) <- attribute_list[[i]][vector]
        } else {
            attr(result, name_attribute[i]) <- attribute_list[[i]]
        }

    }

    return(result)


}


#' @title Output with listed documents
#' @description Behaves like \code{\link[base]{cat}}, but it first automatically unlists the exam to print the document.
#'
#' Since the document is kept as a tree of lists, it simply abstract the idea of outputting the document. with one document.
#'
#'
#' @param FullDocument Document as structure by \code{\link{StructureDocument}}
#' @param sep The separation character(s) between each line.
#' @param ... all extra arguments get passed along to the command "\code{\link[base]{cat}}"
#'
#'
#' @export
#' @examples
#' catDocument(TexExamRandomizer::testdoc)
#'
#'
catDocument <- function(FullDocument, sep = "\n", ...){
    cat(unlist(FullDocument),sep = sep, ...)
}


# @example
# catDocument(Document, sep = "\n", file = stdout())
#
#





#' @title Apply function within a folder
#'
#' @description It executes the function \code{fun} by first switching directories temporarily to the folder \code{folder} and then returning to the working directory.
#'
#' @param folder The folder of execution that the function is switched to before executing \code{fun}
#' @param fun Function to be executed from the relative path
#' @param ... Options to be passed to \code{fun}
#'
#' @return The return value of \code{fun(...)}
#'
#'
#' @export
#' @examples
#' list.files()
#' fun_from_folder(system.file("data", package = "TexExamRandomizer"), list.files)
#' list.files()
fun_from_folder <- function(folder, fun, ...) {
    owd <- getwd() #Original Working Directory
    if (!is.null(owd)) {
        on.exit(setwd(owd), add = TRUE)
        setwd(folder)
        #After the on.exit, if it causes an error setwd(owd) is still called, do not put the set before the on.exit....
    } else {
        warning("Current working directory is unknown, can't change the compilation directory and return to the same directory afterwards...")
    }

    return(fun(...))


}


