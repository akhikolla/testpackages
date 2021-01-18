
#' @rdname match
#' @export
match <- function(x, table, ...) {
  UseMethod("match")
}

#' @rdname match
#' @export
match.default <- function(x, table, ...) {
  base::match(x, table, ...)
}


#' Value matching
#' 
#' @param x \code{\link{lvec}} of values to be matched
#' @param table vector of values in which to look for matches.
#' @param na_incomparable can NA's and NaN's be matched. 
#' @param ... optional arguments passed to and from other methods.
#' 
#' @return 
#' Returns a numeric \code{\link{lvec}} of the same length as \code{x} with the
#' corresponding indices of records in table with the same value. When no match
#' in \code{table} is found, \code{NA} is returned for the corresponding record.
#'
#' @rdname match
#' @useDynLib ldat
#' @export
match.lvec <- function(x, table, na_incomparable = FALSE, ...) {
  if (!is_lvec(table)) table <- as_lvec(table)
  ox <- order(x)
  otable <- order(table)
  structure(lmatch_cpp(x, ox, table, otable, na_incomparable), 
    class = "lvec")
  #structure(.Call("lmatch", PACKAGE = "ldat", x, ox, table, otable, na_incomparable), 
    #class = "lvec")
}

