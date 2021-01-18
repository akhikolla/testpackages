################################################################################

#' Class mmapchar
#'
#' A reference class for storing and accessing matrix-like data stored on disk 
#' in files containing only characters (digits) separated by a character.
#'
#' @examples
#' test_file <- system.file("testdata/test-windows.txt", package = "mmapcharr")
#' test <- mmapchar(test_file, code = mmapcharr:::CODE_012)
#' test[, 1:3]
#' test[]
#' readLines(test_file)
#' 
#' @exportClass mmapchar
#'
mmapchar_RC <- methods::setRefClass(
  
  "mmapchar",
  
  fields = list(
    extptr = "externalptr",
    nrow = "integer",
    ncol = "integer",
    nextra = "integer",
    backingfile = "character",
    code = "vector",
    
    address = function() {
      if (identical(.self$extptr, new("externalptr"))) { # nil
        .self$extptr <- charSepXPtr(.self$backingfile,
                                    .self$nrow,
                                    .self$ncol,
                                    .self$nextra)
      }
      .self$extptr
    }
  ),
  
  methods = list(
    initialize = function(file, d, code) {
      
      if (length(code) != 256) stop2("Parameter 'code' must be of length 256.")
      
      .self$backingfile <- file
      .self$nrow        <- d[["nrow"]]
      .self$ncol        <- d[["ncol"]]
      .self$nextra      <- d[["nextra"]]
      .self$code        <- code
      
      .self$address  # connect once
    },
    
    dim = function() {
      c(nrow = .self$nrow, ncol = .self$ncol, nextra = .self$nextra)
    },
    
    copy = function(code = .self$code) {
      new(Class = "mmapchar", file = .self$backingfile, d = .self$dim(), 
          code = code)
    },
    
    show = function() {
      cat(sprintf("A mmapchar with %d rows and %d columns.",
                  .self$nrow, .self$ncol))
      invisible(.self)
    }
  )
)
mmapchar_RC$lock("backingfile", "nrow", "ncol", "nextra")

################################################################################

#' Wrapper constructor for class `mmapchar`.
#'
#' @param file Path of the file.
#' @param code Integer vector of size 256 to access integers instead of
#'   `rawToChar(as.raw(0:255), multiple = TRUE)`. 
#'   See `mmapcharr:::CODE_012` and `mmapcharr:::CODE_DIGITS`.
#'
#' @rdname mmapchar-class
#'
#' @export
#'
mmapchar <- function(file, code) {
  new(Class = "mmapchar", file = file, d = dim_file(file), code = code)
}

################################################################################

#' Methods for the mmapchar class
#'
#' @name mmapchar-methods
#'
#' @rdname mmapchar-methods
NULL

#' Accessor methods for class `mmapchar`. You can use positive and negative indices,
#' logical indices (that are recycled) and also a matrix of indices (but only
#' positive ones).
#'
#' @param x A [mmapchar][mmapchar-class] object.
#' @param i A vector of indices (or nothing). You can use positive and negative
#'   indices, logical indices (that are recycled) and also a matrix of indices
#'   (but only positive ones).
#' @param j A vector of indices (or nothing). You can use positive and negative
#'   indices, logical indices (that are recycled).
#' @param ... Not used. Just to make [nargs] works.
#' @param drop Whether to delete the dimensions of a matrix which have
#'   one dimension equals to 1.
#'
#' @rdname mmapchar-methods
#'
#' @include extract.R
#'
#' @export
#'
setMethod(
  '[', signature(x = "mmapchar"),
  Extract(extract_vector = extractVec, extract_matrix = extractMat)
)

################################################################################

#' Dimension and type methods for class `mmapchar`.
#'
#' @rdname mmapchar-methods
#' @export
setMethod("dim",    signature(x = "mmapchar"), function(x) c(x$nrow, x$ncol))

#' @rdname mmapchar-methods
#' @export
setMethod("length", signature(x = "mmapchar"), function(x) prod(dim(x)))

################################################################################