#' @useDynLib gRim

## Vanilla R imports (and exports)
## -------------------------------

#' @importFrom stats addmargins as.formula cov.wt fitted formula
#'     ftable getCall logLik loglin na.omit pchisq pf pnorm r2dtable
#'     terms update update.formula xtabs
#' @importFrom utils combn str
#' 
#' @importMethodsFrom stats4 plot
#' 

## To make available in vignette 
#' @importFrom magrittr   "%>%"
#' @export "%>%" 

## Miscellaneous
## -------------

#' @importFrom Rcpp evalCpp
#'
#' @import methods
#' @import gRbase
#'
#' @importFrom gRain propagateLS
#' @importFrom igraph igraph.to.graphNEL igraph.from.graphNEL
#'     get.adjacency V "V<-" E "E<-" is.directed layout.lgl
#'     plot.igraph graph.adjacency is.dag

## Bioconductor imports/exports
## ----------------------------

#' @importClassesFrom graph graphNEL
#' @importFrom graph edges nodes

#' @importMethodsFrom Rgraphviz plot

.dumfunction_afterimportFrom <- function(){}

