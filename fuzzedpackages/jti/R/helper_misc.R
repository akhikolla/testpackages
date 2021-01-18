## ---------------------------------------------------------
##                NON-EXPORTED HELPERS
## ---------------------------------------------------------

## MAPS
.map_chr     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = character(1), ...)
.map_int     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = integer(1), ...)
.map_dbl     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = numeric(1), ...)
.map_lgl     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = logical(1), ...)
.map_lst     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = list(), ...)

## STRINGS
## str_rem <- function(s, pos) {
##   # Vectorized removal of substrings
##   # s: character vector
##   .map_chr(strsplit(s, ""), function(x) {
##     paste0(x[-pos], collapse = "")
##   })
## }

## GRAPHS
as_adj_lst <- function(A) {
  names_ <- colnames(A)
  out <- lapply(seq_along(names_), function(r) {
    names_[as.logical(A[, r])]
  })
  names(out) <- names_
  out
}

as_adj_mat <- function(adj) {
  # TODO: Convert to c++ function using arma::Mat
  stopifnot(length(names(adj)) == length(adj))
  names_ <- names(adj)
  N      <- length(names_)
  A      <- matrix(0L, nrow = N, ncol = N, dimnames = list(names_, names_))
  for (d in seq_along(names_)) {
    idx <- match(adj[[d]], names_)
    A[idx, d] <- 1L
    A[d, idx] <- 1L
  }
  A
}

## MISC
push <- function(l, el, name = NULL) {
  # TODO: if el is a named list, we must take this into account
  c(l, structure(list(el), names = name))
}
