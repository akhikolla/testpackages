## This is a dummy function to let roxygen2 @import the tm, igraph and Matrix package 
## This way, roxygen2 will add them to namespace

#' @import tm 
#' @import Matrix
#' @import igraph
#' @importFrom Rcpp evalCpp
#' @useDynLib RNewsflow
importglobals <- function(x) NULL

