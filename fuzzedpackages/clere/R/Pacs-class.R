#' \code{\linkS4class{Pacs}} class
#' 
#' This class contains all the input parameters to run CLERE.
#' 
#' \describe{ 
#'   \item{Y}{[numeric]: The vector of observed responses - size \code{n}.} 
#'   \item{X}{[matrix]: The matrix of predictors - size \code{n} rows and \code{p} columns.} 
#'   \item{lambda}{[numeric]: A non-negative penalty term that controls simultaneouly clusetering and sparsity.}
#'   \item{betaInput}{[numeric]: A vector of initial guess of the model parameters. The authors suggest to use coefficients obtained after fitting a ridge regression with the shrinkage parameter selected using AIC criterion.}
#'   \item{epsPACS}{[numeric]: A tolerance threshold that control the convergence of the algroithm. The default value fixed in Bondell's initial script is 1e-5.}   \item{nItMax}{[integer]: Maximum number of iterations in the algorithm.} 
#'   \item{a0}{[numeric]: Fitted intercept.} \item{K}{[integer]: Model dimensionality.} 
#' }
#' 
#' @name Pacs-class
#' @docType class
#' @section Methods: 
#' \describe{ 
#'   \item{object["slotName"]:}{Get the value of the field \code{slotName}.} 
#'   \item{object["slotName"]<-value:}{Set \code{value} to the field \code{slotName}.}
#' }
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}}, \code{\linkS4class{Pacs}} \cr 
#' Methods : \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}}, \code{\link{fitPacs}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}, \code{\link{algoComp}}
#' 
#' @keywords Pacs Clere class methods method
#' @export
#' 
methods::setClass(
  Class = "Pacs",
  representation = methods::representation(
    y = "numeric",
    x = "matrix",
    n = "integer",
    p = "integer",
    nItMax = "integer",
    lambda = "numeric",
    epsPACS = "numeric",
    betaInput = "numeric",
    betaOutput = "numeric",
    a0 = "numeric",
    K = "integer"
  ),
  prototype = methods::prototype(
    y = numeric(),
    x = matrix(),
    n = integer(),
    p = integer(),
    nItMax = integer(),
    lambda = numeric(),
    epsPACS = numeric(),
    betaInput = numeric(),
    betaOutput = numeric(),
    a0 = numeric(),
    K = integer()
  )
)
