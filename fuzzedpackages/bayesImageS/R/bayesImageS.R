#' Package bayesImageS
#'
#' Bayesian methods for segmentation of 2D and 3D images, such as computed
#' tomography and satellite remote sensing. This package implements image
#' analysis using the hidden Potts model with external field prior. Latent labels
#' are sampled using chequerboard updating or Swendsen-Wang. Algorithms for the 
#' smoothing parameter include pseudolikelihood, path sampling, the exchange
#' algorithm, and approximate Bayesian computation (ABC-MCMC and ABC-SMC).
#' 
#' @author
#' M. T. Moores and K. Mengersen
#' with additional code contributed by D. Feng
#' 
#' Maintainer: Matt Moores <M.T.Moores@@warwick.ac.uk>
#' 
#' @references
#' Moores, M. T.; Pettitt, A. N. & Mengersen, K. (2015) "Scalable Bayesian inference for the inverse temperature of a hidden Potts model" \href{http://arxiv.org/abs/1503.08066}{arXiv:1503.08066}
#' 
#' Moores, M. T.; Drovandi, C. C.; Mengersen, K. & Robert, C. P. (2015) "Pre-processing for approximate Bayesian computation in image analysis" \emph{Statistics & Computing} \bold{25}(1), 23--33, DOI: \href{http://dx.doi.org/10.1007/s11222-014-9525-6}{10.1007/s11222-014-9525-6}
#' 
#' Moores, M. T.; Hargrave, C. E.; Deegan, T.; Poulsen, M.; Harden, F. & Mengersen, K. (2015) "An external field prior for the hidden Potts model, with application to cone-beam computed tomography" \emph{Computational Statistics & Data Analysis} \bold{86}, 27--41, DOI: \href{http://dx.doi.org/10.1016/j.csda.2014.12.001}{10.1016/j.csda.2014.12.001}
#' 
#' Feng, D. (2008) "Bayesian Hidden Markov Normal Mixture Models with Application to MRI Tissue Classification" \emph{Ph. D. Dissertation, The University of Iowa} 
#'
#' @seealso
#' \code{vignette(package="bayesImageS")}
#'
#' @useDynLib bayesImageS
#' @importFrom Rcpp sourceCpp
#' @docType package
#' @name bayesImageS
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("bayesImageS", libpath)
}