#' @useDynLib CASMAP
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom methods new
NULL

#' Global variables environment
#'
#' An environment to store a few global variables. Internal.
#'
#' @keywords internal
CASMAPenv <- new.env(parent=emptyenv())
assign("regionGWASString", "regionGWAS", envir=CASMAPenv)
assign("higherOrderEpistasisString", "higherOrderEpistasis", envir=CASMAPenv)
assign("minModeLength", 3, envir=CASMAPenv)
