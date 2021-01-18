
#' chunker initializer
#' @aliases initialize,chunker-method
#' @keywords internal

setMethod( "initialize", "chunker", 
           function(.Object, path_, sep_, quoted_, has_colnames_, 
                    has_rownames_, chunksize_, data_format_,  columns_classes_) {
             
             path_ <- normalizePath(path_)
             
             if(data_format_ == "matrix") {
               .Object@pointer  <- chunker__new_matrix(path_, 
                                                      sep_, 
                                                      quoted_,
                                                      has_colnames_, 
                                                      has_rownames_, 
                                                      chunksize_)
             } else {
               .Object@pointer  <- chunker__new_data_frame(path_, 
                                                           sep_,
                                                           quoted_,
                                                           has_colnames_, 
                                                           has_rownames_,
                                                           chunksize_, 
                                                           columns_classes_)
             }
             .Object
           })

#' Manipulation methods for chunker objects
#' @name chunker methods
#' @description chunker objects can be manipulated with the following methods:
#' \enumerate{
#' \item{\bold{next_chunk}}{: allows to read the next chunk of a chunker object}
#' \item{\bold{get_table}}{: retrieve the current data chunk contained in the object}
#' }
#' 
#' In addition, the following information can be retrieved from chunker objects:
#' \enumerate{
#' \item{\bold{get_completed}}{: get the number of rows already read}
#' \item{\bold{get_colnames}}{: get column names of the chunker object}
#' }
#' 
#' @details See \code{\link{chunker}} for examples.
#' 
#' @param obj object of class chunker
#' @rdname chunker-methods
NULL


#' next_chunk
#' @name next_chunk
#' @description NULL
#' @rdname chunker-methods
#' @export
setGeneric("next_chunk", function(obj) standardGeneric("next_chunk"))


#' @rdname chunker-methods
#' @aliases next_chunk,chunker-methods

setMethod("next_chunk", "chunker", function(obj) {
  chunker__next_chunk(obj@pointer)
})


#' get_table
#' @name get_table
#' @description NULL
#' @rdname chunker-methods
#' @export
setGeneric("get_table", function(obj) standardGeneric("get_table"))


#' @aliases get_dataframe, chunker-methods
#' @rdname chunker-methods

setMethod("get_table", "chunker", function(obj) {
  what_is <- get_type(obj) 
  if(what_is == "data.frame") {
    chunker__get_dataframe(obj@pointer)
  } else {
    chunker__get_matrix(obj@pointer)
  }
})


#' get_colnames
#' @name get_colnames
#' @rdname chunker-methods
#' @description NULL
#' @export

setGeneric("get_colnames", function(obj) standardGeneric("get_colnames"))


#' @aliases get_colnames,chunker-method
#' @rdname chunker-methods

setMethod("get_colnames", "chunker", function(obj) {
  chunker__get_colnames(obj@pointer)
})


#' get_completed
#' @name get_completed
#' @description NULL
#' @rdname chunker-methods
#' @export

setGeneric("get_completed", function(obj) standardGeneric("get_completed"))


#' @aliases get_completed,chunker-methods
#' @rdname chunker-methods

setMethod("get_completed", "chunker", function(obj) {
  chunker__get_completed(obj@pointer)
})

#' get_total
#' @name get_total
#' @description NULL
#' @rdname chunker-methods

setGeneric("get_total", function(obj) standardGeneric("get_total"))


#' @aliases get_total,chunker-methods
#' @rdname chunker-methods

setMethod("get_total", "chunker", function(obj) {
  chunker__get_total(obj@pointer)
})


#' get_type
#' @name get_type
#' @description NULL
#' @rdname chunker-methods
#' @export

setGeneric("get_type", function(obj) standardGeneric("get_type"))


#' @aliases get_type,chunker-methods
#' @rdname chunker-methods

setMethod("get_type", "chunker", function(obj) {
  chunker__get_type(obj@pointer)
})

#' get_attr
#' @name get_attr
#' @description NULL
#' @rdname chunker-methods
#' @export

setGeneric("get_attr", function(obj) standardGeneric("get_attr"))


#' @aliases get_attr,chunker-methods
#' @rdname chunker-methods

setMethod("get_attr", "chunker", function(obj) {
  obj@attr
})

#' @aliases print,chunker-methods
#' @keywords internal

setMethod("show", "chunker", function(object) {
  total_lines <- get_total(object)
  completed_lines <- get_completed(object)
  cat("\n>---- chunker object ---->\n\n")
  cat("Path: ", object@attr$path, "\n")
  cat("Chunk size: ", object@attr$chunksize, "| Sep: ", deparse(object@attr$sep), "| Quoted: ", object@attr$quoted, "\n")
  cat("Colnames: ", object@attr$has_colnames, "| Rownames: ", object@attr$has_rownames, "\n")
  cat("Data imported as: ", object@attr$data_format, "\n")
  cat("Total lines in file: ", ifelse (object@attr$data_format == "data.frame", total_lines, "unknown"), "\n")
  cat("Lines completed: ", completed_lines)
  cat(" ", ifelse (object@attr$data_format == "data.frame", paste0("(", round(100 * completed_lines/total_lines), "%)\n"), "\n"))
  cat("Additional information: use the function 'get_attr' with this object", "\n\n")
})


#-------------
# Add documentation for matrix2df

#' matrix2df
#' @description conversion from matrix to data frame
#' @param x A matrix
#' @export
NULL
