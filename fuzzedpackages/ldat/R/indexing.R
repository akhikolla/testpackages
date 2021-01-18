
#' Indexing of lvec objects
#'
#' @param x an object of type \code{\link{lvec}}
#' @param i an index vector. See \code{\link{lget}}.
#' @param j a selection of columns (a character, numeric or logical vector).
#' @param range an range of indices. See \code{\link{lget}}.
#' @param value new values. See \code{\link{lget}}.
#' @param clone \code{\link{clone}} columns when selecting only columns. 
#' @param drop ignored; included for compatability with \code{data.frame}.
#'
#' @details
#' These functions are a wrapper around \code{\link{lget}} and 
#' \code{\link{lset}}.
#'
#' @rdname indexing
#' @import lvec
#' @export
`[.lvec` <- function(x, i = NULL, range = NULL) {
  lget(x, index = i, range = range)
}

#' @rdname indexing
#' @import lvec
#' @export
`[<-.lvec` <- function(x, i, range, value) {
  if (!missing(range)) {
    lset(x, range = range, values = value)
  } else {
    lset(x, index = i, values = value)
  }
}


#' @rdname indexing
#' @export
`[.ldat` <- function(x, i, j, drop = FALSE, range = NULL, clone = TRUE) {
  nindices <- nargs() - 1L
  if (!missing(clone)) nindices <- nindices - 1L
  if (!missing(drop)) nindices <- nindices - 1L
  if (!missing(drop) && drop)
    warning("'drop = TRUE' argument will be ignored.")

  if (nindices == 0) {
    res <- if (clone) clone(x) else x
  } else if (nindices == 1 && missing(range)) {
    res <- unclass(x)[i]
    if (clone) res <- lapply(res, clone)
    structure(res, class = "ldat")
  } else if (nindices == 1 && !missing(range)) {
    lget(x, range = range)
  } else if (nindices == 2 && !missing(range)) {
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- i
    x <- x[j, clone = FALSE]
    lget(x, range = range)
  } else if (nindices == 2 && missing(range)) {
    if (missing(i)) i <- TRUE
    if (!missing(j)) x <- x[j, clone = FALSE]
    lget(x, index = i)
  } else {
    stop("Invalid indices.")
  }
}

#' @rdname indexing
#' @export
`[<-.ldat` <- function(x, i, range, value) {
  if (!missing(range)) {
    lset(x, range = range, values = value)
  } else {
    lset(x, index = i, values = value)
  }
}
