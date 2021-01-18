


#' @export
dim.ldat <- function(x) {
  ncol <- length(x)
  nrow <- if (ncol == 0) 0 else length(x[[1]])
  c(nrow, ncol)
}
