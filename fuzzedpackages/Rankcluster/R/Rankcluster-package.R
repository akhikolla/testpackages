#' @import Rcpp methods
#' @useDynLib Rankcluster
#'
#' @title Model-Based Clustering for Multivariate Partial Ranking Data
#' @docType package
#' @aliases Rankcluster-package
#' @name Rankcluster-package
#' 
#' @description This package proposes a model-based clustering algorithm for ranking data.
#' Multivariate rankings as well as partial rankings are taken into account.
#' This algorithm is based on an extension of the Insertion Sorting Rank (ISR) model for ranking data, which is a meaningful
#' and effective model parametrized by a position parameter (the modal ranking, quoted by mu) and a dispersion parameter (quoted by pi).
#' The heterogeneity of the rank population is modelled by a mixture of ISR, whereas conditional independence assumption is considered for multivariate rankings.
#'
#' @details
#' The main function is \link{rankclust}.
#' See vignettes for detailled examples: \code{RShowDoc("dataFormat", package = "Rankcluster")} and \code{RShowDoc("Rankcluster", package = "Rankcluster")}
#' 
#'
#' @references   [1] C.Biernacki and J.Jacques (2013), A generative model for rank data based on sorting algorithm, Computational Statistics and Data Analysis, 58, 162-176.
#'
#' [2] J.Jacques and C.Biernacki (2012), Model-based clustering for multivariate partial ranking data, Inria Research Report n 8113.
#'
#' @author Maintainer: Quentin Grimonprez <quentin.grimonprez@@inria.fr>
#'
#' @examples
#' # see vignettes
#' # RShowDoc("dataFormat", package = "Rankcluster")
#' # RShowDoc("Rankcluster", package = "Rankcluster")
#'
#' # main function of the package for run the algorithm
#' data(big4)
#' result <- rankclust(big4$data, K = 2, m = big4$m, Ql = 200, Bl = 100, maxTry = 2)
#' 
#' if(result@@convergence)
#' { 
#'   summary(result)
#' 
#'   partition <- result[2]@@partition
#'   tik <- result[2]@@tik
#' }
#' 
#' @keywords package
NULL
