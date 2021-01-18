## SFS.R -- R-level functionality for SFS
## Currently only used to let roxygen (and devtools::check) generate the
## NAMESPACE file
#' @useDynLib SFS, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL
