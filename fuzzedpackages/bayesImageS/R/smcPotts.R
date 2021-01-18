#' Fit the hidden Potts model using approximate Bayesian computation with sequential Monte Carlo (ABC-SMC).
#' 
#' @param y A vector of observed pixel data.
#' @param neighbors A matrix of all neighbors in the lattice, one row per pixel.
#' @param blocks A list of pixel indices, dividing the lattice into independent blocks.
#' @param param A list of options for the ABC-SMC algorithm.
#' @param priors A list of priors for the parameters of the model.
#' @return A matrix containing SMC samples for the parameters of the Potts model.
#' @export
smcPotts <- function(y, neighbors, blocks, param=list(npart=10000, nstat=50), priors=NULL) {
  result <- .Call( "smcPotts", y, neighbors, blocks, param, priors, PACKAGE = "bayesImageS")
}

#' Initialize the ABC algorithm using the method of Sedki et al. (2013)
#' 
#' @param y A vector of observed pixel data.
#' @param neighbors A matrix of all neighbours in the lattice, one row per pixel.
#' @param blocks A list of pixel indices, dividing the lattice into independent blocks.
#' @param param A list of options for the ABC-SMC algorithm.
#' @param priors A list of priors for the parameters of the model.
#' @return A matrix containing SMC samples for the parameters of the Potts model.
#' @export
#' @references
#' Sedki, M.; Pudlo, P.; Marin, J.-M.; Robert, C. P. & Cornuet, J.-M. (2013) "Efficient learning in ABC algorithms" \href{http://arxiv.org/abs/1210.1388}{arXiv:1210.1388}
initSedki <- function(y, neighbors, blocks, param=list(npart=10000), priors=NULL) {
  result <- .Call( "initSedki", y, neighbors, blocks, param, priors, PACKAGE = "bayesImageS")
}

#' Test the residual resampling algorithm.
#' 
#' @param values A vector of SMC particles.
#' @param weights A vector of importance weights for each particle.
#' @param pseudo A matrix of pseudo-data for each particle.
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{beta}}{A vector of resampled particles.}
#'   \item{\code{wt}}{The new importance weights, after resampling.}
#'   \item{\code{pseudo}}{A matrix of pseudo-data for each particle.}
#'   \item{\code{idx}}{The indices of the parents of the resampled particles.}
#'   }
#' @export
#' @references 
#'  Liu, J. S. & Chen, R. (1998) "Sequential Monte Carlo Methods for Dynamic Systems"
#'  \emph{J. Am. Stat. Assoc.} \bold{93}(443): 1032--1044, DOI: \href{https://doi.org/10.1080/01621459.1998.10473765}{10.1080/01621459.1998.10473765}
testResample <- function(values, weights, pseudo) {
  result <- .Call( "testResample", values, weights, pseudo, PACKAGE = "bayesImageS")
}

#' Calculate the distribution of the Potts model using a brute force algorithm.
#' 
#' \bold{Warning}: this algorithm is O\eqn{(k^n)} and therefore will not scale for
#' \eqn{k^n > 2^{31} - 1}
#' 
#' @param neighbors A matrix of all neighbours in the lattice, one row per pixel.
#' @param blocks A list of pixel indices, dividing the lattice into independent blocks.
#' @param k The number of unique labels.
#' @param beta The inverse temperature parameter of the Potts model.
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{expectation}}{The exact mean of the sufficient statistic.}
#'   \item{\code{variance}}{The exact variance of the sufficient statistic.}
#'   \item{\code{exp_PL}}{Pseudo-likelihood (PL) approximation of the expectation of S(z).}
#'   \item{\code{var_PL}}{PL approx. of the variance of the sufficient statistic.}
#'   }
#' @export
exactPotts <- function(neighbors, blocks, k, beta) {
  result <- .Call( "exactPotts", neighbors, blocks, k, beta, PACKAGE = "bayesImageS")
}
