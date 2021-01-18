#' Functions and C++ Header Files for Computation on Manifolds
#' 
#' We provide a number of algorithms to estimate fundamental statistics including 
#' Fréchet mean and geometric median for manifold-valued data. Also, C++ header files are contained that implement elementary operations 
#' on manifolds such as Sphere, Grassmann, and others. See Bhattacharya and Bhattacharya (2012) <doi:10.1017/CBO9781139094764> 
#' if you are interested in statistics on manifolds, and Absil et al (2007) <isbn:978-0-691-13298-3> on computational 
#' aspects of optimization on matrix manifolds.
#' 
#' @docType package
#' @name RiemBase-package
#' @import Rdpack
#' @importFrom pracma flipud
#' @importFrom stats quantile rnorm
#' @importFrom parallel detectCores
#' @importFrom utils getFromNamespace packageVersion
#' @importFrom Rcpp evalCpp
#' @useDynLib RiemBase
NULL
