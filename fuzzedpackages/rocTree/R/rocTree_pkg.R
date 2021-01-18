 #' rocTree:Receiver Operating Characteristic (ROC)-Guided Classification Survival Tree and Ensemble.
#'
#' The \code{rocTree} package uses a Receiver Operating Characteristic (ROC) guided classification
#' algorithm to grow prune survival trees and ensemble.
#' 
#'
#' @aliases rocTree-package
#' @section Introduction:
#' The \code{rocTree} package provides implementations to a unified framework for
#' tree-structured analysis with censored survival outcomes.
#' Different from many existing tree building algorithms,
#' the \code{rocTree} package incorporates time-dependent covariates by constructing
#' a time-invariant partition scheme on the survivor population.
#' The partition-based risk prediction function is constructed using an algorithm guided by
#' the Receiver Operating Characteristic (ROC) curve.
#' The generalized time-dependent ROC curves for survival trees show that the
#' target hazard function yields the highest ROC curve.
#' The optimality of the target hazard function motivates us to use a weighted average of the
#' time-dependent area under the curve on a set of time points to evaluate the prediction
#' performance of survival trees and to guide splitting and pruning.
#' Moreover, the \code{rocTree} package also offers a novel ensemble algorithm,
#' where the ensemble is on unbiased martingale estimating equations.
#' 
#' @section Methods:
#' The package contains functions to construct ROC-guided survival trees and ensemble through
#' the main function \code{\link{rocTree}}.
#' 
#' @seealso \code{\link{rocTree}}
#' @docType package
#' 
#' @importFrom stats model.extract model.matrix model.response
#' @importFrom utils tail
#' @importFrom survival survfit Surv
#' @importFrom parallel detectCores makeCluster setDefaultCluster clusterExport stopCluster
#' @importFrom parallel parSapply parLapply
#' @importFrom graphics legend lines plot
#' @importFrom data.tree Node ToDataFrameTree ToDiagrammeRGraph SetGraphStyle SetNodeStyle
#' @importFrom DiagrammeR render_graph %>% export_graph
#' @importFrom Rcpp sourceCpp
#' 
#' @useDynLib rocTree 
"_PACKAGE"
NULL


## @useDynLib rocTree, .registration = TRUE
