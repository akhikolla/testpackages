#' @useDynLib sprintr
#' @importFrom Rcpp sourceCpp
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("sprintr", libpath)
}
