## ---------------------------------------------------------
##                NON-EXPORTED HELPERS
## ---------------------------------------------------------
## MAPS
.map_dbl     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = numeric(1), ...)
.map_lgl     <- function(x, fun, ...) vapply(X = x, FUN = fun, FUN.VALUE = logical(1), ...)

## STRINGS
.split_chars <- function(x) unlist(strsplit(x, ""))

## MISC
only_single_chars <- function(A) {
  for (i in seq_along(nrow(A))) {
    for (j in seq_along(ncol(A)))
      if ( nchar(A[i,j]) != 1L ) return(FALSE)
  }
  return(TRUE)
}


## ---------------------------------------------------------
##                     EXPORTED HELPERS
## ---------------------------------------------------------
#' Convert discrete values into a single character representation
#'
#' Convert all values in a data frame or matrix of characters to a single
#' character representation
#'
#' @param x Data frame or matrix of characters
#' @examples
#' d <- data.frame(x = c("11", "2"), y = c("2", "11"))
#' to_single_chars(d)
#' @export
to_single_chars <- function(x) {
  # Implicitly assumes that no columns has more than length(chars) = 62 unique levels
  # Consider saving the olde levels so we can retrive them again easily later
  apply(x, 2, function(z) {
    f <- as.factor(z)
    chars <- c(letters, LETTERS, 0:9)
    levels(f) <- chars[1:length(levels(f))]
    as.character(f)
  })
}
