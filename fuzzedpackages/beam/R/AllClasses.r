################
## All Classes
################

#' @name beam-class
#' @aliases beam-class
#'
#' @title Class beam
#'
#' @description An S4 class representing the output of the \code{\link{beam}} function.
#'
#' @slot table matrix. A data.frame containing marginal and/or partial correlation estimates, Bayes factors and tail probabilities for each edge.
#' @slot deltaOpt numeric. Empirical Bayes estimate of hyperpaprameter delta.
#' @slot alphaOpt numeric. Empirical Bayes estimate of hyperpaprameter alpha.
#' @slot dimX numeric. Dimension of the input data matrix X.
#' @slot type character. Input argument.
#' @slot varlabs character. Column labels of X.
#' @slot gridAlpha matrix. A matrix containing the log-marginal likelihood of the Gaussian conjugate model as a function of a grid of values of alpha and delta.
#' @slot valOpt numeric. Maximum value of the log-marginal likelihood of the Gaussian conjugate model.
#' @slot return.only character. Input argument.
#' @slot time numeric. Running time (in seconds).
#' @slot TinvStdev numeric. Square root of partial variances.
#' @slot s numeric. Sample variances.
#'
#' @method
#'  \item{print}{Print the object information}
#'  \item{show}{Print the object information}
#'  \item{summary}{Print summary of the object information.}
#'  \item{marg}{Return a data.frame containing marginal correlation estimates, Bayes factors and tail probabilities for all edges.}
#'  \item{cond}{Return a data.frame containing partial correlation estimates, Bayes factors and tail probabilities for all edges.}
#'  \item{cormat}{Return the marginal correlation matrix.}
#'  \item{pcormat}{Return the partial correlation matrix.}
#'  \item{plotML}{Plot the log-marginal likelihood of the Gaussian conjugate model as a function of shrinkage parameter alpha.}
#'  \item{plotHeatmap}{Plot heatmap of marginal (upper triangle) and/or partial (lower triangle) correlation estimates.}
#'  \item{bgraph}{Extract the bidirected graph (or marginal inddependence graph) as an igraph object.}
#'  \item{ugraph}{Extract the bidirected graph (or marginal inddependence graph) as an igraph object.}
#'
#' @author Gwenael G.R. Leday and Ilaria Speranza
#'
setClass("beam",
         representation(table = "matrix",
                        deltaOpt = "numeric",
                        alphaOpt = "numeric",
                        dimX = "numeric",
                        type = "character",
                        varlabs = "character",
                        gridAlpha = "matrix",
                        valOpt = "numeric",
                        return.only = "character",
                        time = "numeric",
                        TinvStdev = "numeric",
                        s = "numeric")
)

#' @name beam.select-class
#' @aliases beam.select-class
#'
#' @title Class beam.select
#'
#' @description An S4 class representing the output of the \code{\link{beam.select}} function.
#'
#' @slot marginal data.frame. A data.frame containing the marginal correlation estimates, Bayes factors and tail probabilities for the selected edges only.
#' @slot conditional data.frame. A data.frame containing the partial correlation estimates, Bayes factors and tail probabilities for the selected edges only.
#' @slot dimX numeric. Dimension of the imput data matrix X.
#' @slot type character. Input type (marginal, conditional or both)
#' @slot varlabs character. Column labels of X.
#' @slot alphaOpt numeric. Empirical Bayes estimates of hyperpaprameter alpha.
#' @slot gridAlpha matrix. A matrix containing the log-marginal likelihood of the Gaussian conjugate model as a function of a grid of values of alpha and delta.
#' @slot valOpt numeric. Maximum value of the log-marginal likelihood of the Gaussian conjugate model
#' @slot method character. Input method.
#' @slot thres numeric. Input threshold

#' @method
#'  \item{print}{Print the object information}
#'  \item{show}{Print the object information}
#'  \item{summary}{Print summary of the object information.}
#'  \item{marg}{Return a data.frame containing marginal correlation estimates, Bayes factors and tail probabilities for selected edges.}
#'  \item{cond}{Return a data.frame containing partial correlation estimates, Bayes factors and tail probabilities for selected edges.}
#'  \item{cormat}{Return the (possibly sparse) marginal correlation matrix.}
#'  \item{pcormat}{Return the (possibly sparse) partial correlation matrix.}
#'  \item{plotML}{Plot the log-marginal likelihood of the Gaussian conjugate model as a function of shrinkage parameter alpha.}
#'  \item{plotHeatmap}{Plot heatmap of marginal (upper triangle) and/or partial (lower triangle) correlation estimates.}
#'  \item{bgraph}{Extract the bidirected graph (or marginal inddependence graph) as an igraph object.}
#'  \item{ugraph}{Extract the bidirected graph (or marginal inddependence graph) as an igraph object.}
#'
#' @author Gwenael G.R. Leday and Ilaria Speranza
#'
setClass("beam.select",
         representation(marginal = "data.frame",
                        conditional = "data.frame",
                        dimX = "numeric",
                        type = "character",
                        varlabs = "character",
                        alphaOpt = "numeric",
                        gridAlpha = "matrix",
                        valOpt = "numeric",
                        method = "character",
                        thres= "numeric"
                        )
)
