
## TODO: tests

#' Create an ldat object
#'
#' This function creates an \code{ldat} object, which behaves similar to a 
#' \code{data.frame} except that its columns are \code{\link{lvec}}. This 
#' allows and \code{ldat} to have an arbitrary large number of rows without
#' running into memory problems. 
#' 
#' @param ... these arguments are of either the form `tag = value' or `value'. 
#'   Each argument becomes a column in the \code{ldat}. All columns are 
#'   required to have the same length. 
#'
#' @details
#' Each of the arguments of \code{ldat} is converted to an \code{\link{lvec}} 
#' when it isn't already and \code{\link{lvec}} using calls to 
#' \code{\link{as_lvec}}. The arguments are required to all have the same 
#' length (unlike \code{\link{data.frame}}).
#' 
#' @return
#' An object of type \code{ldat}. This object is basically a list with 
#' \code{\link{lvec}} objects.
#'
#' @examples
#' # Create ldat object from r-objects
#' a <- ldat(id = 1:20, x = letters[1:20], y = rnorm(20))
#' # this is identical to
#' a <- ldat(id = as_lvec(1:20), x = as_lvec(letters[1:20]), 
#'     y = as_lvec(rnorm(20)))
#'
#' @export
ldat <- function(...) {
  x <- list(...)
  n <- length(x)
  if (n < 1L) {
    return (NULL);
  }
  # Create the column names of the resulting ldat
  vnames <- names(x)
  if (length(vnames) != n) vnames <- character(n)
  object <- as.list(substitute(list(...)))[-1L]
  vnames_deparse <- sapply(object, function(o) deparse(o)[1L])
  vnames  <- ifelse(nzchar(vnames), vnames, vnames_deparse)
  noname <- !nzchar(vnames)
  if (any(noname)) 
    vnames[noname] <- paste("Var", seq_along(vnames), sep = ".")[noname]
  vnames <- make.names(vnames, unique = TRUE)
  # Check each of the columns; if necessary convert to lvec
  len <- NULL
  for (i in seq_len(n)) {
    xi <- as_lvec(x[[i]])
    # Check length; all columns should have same length
    if (is.null(len)) len <- length(xi)
    if (length(xi) != len) {
      # I do not want to change the length of the columns as lvecs are passed
      # by reference and I do not want to change the original value. If the 
      # user wants this he has to do that himself. So generate an error.
      stop("Column '", vnames[i], "' does not have the same length as", 
        " previous columns.")
    }
    x[[i]] <- xi
  }
  names(x) <- vnames
  attr(x, "class") <- "ldat"
  x
}

#' @param x object for which to check if it is of type ldat
#' @rdname ldat
#' @export
is_ldat <- function(x) {
  methods::is(x, 'ldat')
}

