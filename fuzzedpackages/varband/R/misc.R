#' @useDynLib varband
#' @import Rcpp
# #' @import RcppArmadillo
# importing RcppArmadillo here would cause the following note:
# Package in Depends/Imports which should probably only be in LinkingTo: 'RcppArmadillo'
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("varband", libpath)
}
