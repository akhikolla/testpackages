#' ess: Eficient Stepwise Selection in Decomposable Models
#'
#' The class of graphical models is a family of probability
#' distributions for which conditional dependencies can be read off
#' from a graph. If the graph is decomposable, the maximum
#' likelihood estimates of the parameters in the model can be shown
#' to be on exact form. This is what enables ESS to be fast and
#' efficient for model selection in decomposable graphical models.
#' 
#' @importFrom Rcpp sourceCpp
#' @useDynLib ess
"_PACKAGE"
