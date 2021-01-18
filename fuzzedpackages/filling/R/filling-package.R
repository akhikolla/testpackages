#' Matrix Completion, Imputation, and Inpainting Methods
#'
#' Filling in the missing entries of a partially observed data is one of fundamental problems
#' in various disciplines of mathematical science. For many cases, data at our interests
#' have canonical form of matrix in that the problem is posed upon a matrix with missing values
#' to fill in the entries under preset assumptions and models.
#' We provide a collection of methods from multiple disciplines under Matrix Completion, Imputation, and Inpainting.
#' Currently, we have following methods implemented,
#' \tabular{ll}{
#' \emph{Name of Function} \tab \emph{Method} \cr
#' \code{fill.HardImpute} \tab Generalized Spectral Regularization \cr
#' \code{fill.KNNimpute}  \tab Weighted \eqn{K}-nearest Neighbors \cr
#' \code{fill.nuclear} \tab Nuclear Norm Optimization \cr
#' \code{fill.OptSpace} \tab OptSpace \cr
#' \code{fill.simple} \tab Simple Rules of Mean, Median, and Random \cr
#' \code{fill.SoftImpute} \tab Spectral Regularization \cr
#' \code{fill.SVDimpute} \tab Iterative Regression against Right Singular Vectors \cr
#' \code{fill.SVT} \tab Singular Value Thresholding for Nuclear Norm Optimization \cr
#' \code{fill.USVT} \tab Universal Singular Value Thresholding
#' }
#'
#'
#' @docType package
#' @name filling-package
#' @import CVXR
#' @import Rdpack
#' @import nabor
#' @import RSpectra
#' @importFrom utils packageVersion
#' @importFrom ROptSpace OptSpace
#' @importFrom Rcpp evalCpp
#' @importFrom stats median rbinom rnorm lm
#' @useDynLib filling
NULL

# pack <- "filling"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))
