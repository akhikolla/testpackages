#' Matrix Sort
#'
#' Sorts a matrix by first column, breaking ties with second column, breaking those ties with 3rd, etc.
#'
#' @param X matrix to sort.
#' @return Same matrix with rows permuted according to sort order.
#' @keywords internal
colSort <- function(X) {
  if(!is.matrix(X)) return(sort(X))
  X[do.call(order, lapply(data.frame(X), c)),]
}
