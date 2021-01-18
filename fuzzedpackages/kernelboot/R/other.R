
# check if different objects are numeric
# for data.frames and matrix objects check the individual columns

numericColumns <- function(x) UseMethod("numericColumns")

numericColumns.default <- function(x) {
  is.numeric(x)
}

numericColumns.matrix <- function(x) {
  structure(rep(is.numeric(x), ncol(x)), names = colnames(x))
}

numericColumns.data.frame <- function(x) {
  unlist(lapply(x, is.numeric))
}

# check for square matrix

is.square <- function(x) {
  NCOL(x) == NROW(x)
}

# check for diagonal matrix

is.diag <- function(x) {
  dx <- diag(diag(x))
  is.square(x) && all(abs(dx - x) <= .Machine$double.eps)
}

# test for being the vector

is.simple.vector <- function(x) {
  is.atomic(x) && !is.recursive(x) && !is.array(x)
}

# test is all elements are zeros

is.allzeros <- function(x) {
  all(abs(x) <= .Machine$double.eps)
}
