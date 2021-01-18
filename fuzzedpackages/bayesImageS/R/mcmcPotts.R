#' Fit the hidden Potts model using a Markov chain Monte Carlo algorithm.
#' 
#' @param y A vector of observed pixel data.
#' @param neighbors A matrix of all neighbors in the lattice, one row per pixel.
#' @param blocks A list of pixel indices, dividing the lattice into independent blocks.
#' @param priors A list of priors for the parameters of the model.
#' @param mh A list of options for the Metropolis-Hastings algorithm.
#' @param niter The number of iterations of the algorithm to perform.
#' @param nburn The number of iterations to discard as burn-in.
#' @param truth A matrix containing the ground truth for the pixel labels.
#' @return A matrix containing MCMC samples for the parameters of the Potts model.
#' @export
mcmcPotts <- function(y, neighbors, blocks, priors, mh, niter=55000, nburn=5000, truth=NULL) {
  result <- .Call( "mcmcPotts", y, neighbors, blocks, niter, nburn, priors, mh, truth, PACKAGE = "bayesImageS")
}

#' Simulate pixel labels using chequerboard Gibbs sampling.
#' 
#' @param beta The inverse temperature parameter of the Potts model.
#' @param k The number of unique labels.
#' @param neighbors A matrix of all neighbors in the lattice, one row per pixel.
#' @param blocks A list of pixel indices, dividing the lattice into independent blocks.
#' @param niter The number of iterations of the algorithm to perform.
#' @param random Whether to initialize the labels using random or deterministic starting values.
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{alloc}}{An n by k matrix containing the number of times that pixel i was allocated to label j.}
#'   \item{\code{z}}{An \code{(n+1)} by k matrix containing the final sample from the Potts model after niter iterations of chequerboard Gibbs.}
#'   \item{\code{sum}}{An \code{niter} by 1 matrix containing the sum of like neighbors, i.e. the sufficient statistic of the Potts model, at each iteration.}
#'   }
#' @export
#' @examples
#' # Swendsen-Wang for a 2x2 lattice
#' neigh <- matrix(c(5,2,5,3,  1,5,5,4,  5,4,1,5,  3,5,2,5), nrow=4, ncol=4, byrow=TRUE)
#' blocks <- list(c(1,4), c(2,3))
#' res.Gibbs <- mcmcPottsNoData(0.7, 3, neigh, blocks, niter=200)
#' res.Gibbs$z
#' res.Gibbs$sum[200]
mcmcPottsNoData <- function(beta, k, neighbors, blocks, niter=1000, random=TRUE) {
  result <- .Call( "mcmcPottsNoData", beta, k, neighbors, blocks, niter, random, PACKAGE = "bayesImageS")  
}

#' Simulate pixel labels using the Swendsen-Wang algorithm.
#' 
#' The algorithm of Swendsen & Wang (1987) forms clusters of neighbouring pixels,
#' then updates all of the labels within a cluster to the same value. When
#' simulating from the prior, such as a Potts model without an external field,
#' this algorithm is very efficient.
#'
#' @param beta The inverse temperature parameter of the Potts model.
#' @param k The number of unique labels.
#' @param neighbors A matrix of all neighbors in the lattice, one row per pixel.
#' @param blocks A list of pixel indices, dividing the lattice into independent blocks.
#' @param niter The number of iterations of the algorithm to perform.
#' @param random Whether to initialize the labels using random or deterministic starting values.
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{alloc}}{An n by k matrix containing the number of times that pixel i was allocated to label j.}
#'   \item{\code{z}}{An \code{(n+1)} by k matrix containing the final sample from the Potts model after niter iterations of Swendsen-Wang.}
#'   \item{\code{sum}}{An \code{niter} by 1 matrix containing the sum of like neighbors, i.e. the sufficient statistic of the Potts model, at each iteration.}
#'   }
#' @export
#' @references Swendsen, R. H. & Wang, J.-S. (1987) "Nonuniversal critical dynamics in Monte Carlo simulations" \emph{Physical Review Letters} \bold{58}(2), 86--88, DOI: \href{https://doi.org/10.1103/PhysRevLett.58.86}{10.1103/PhysRevLett.58.86}
#' @examples
#' # Swendsen-Wang for a 2x2 lattice
#' neigh <- matrix(c(5,2,5,3,  1,5,5,4,  5,4,1,5,  3,5,2,5), nrow=4, ncol=4, byrow=TRUE)
#' blocks <- list(c(1,4), c(2,3))
#' res.sw <- swNoData(0.7, 3, neigh, blocks, niter=200)
#' res.sw$z
#' res.sw$sum[200]
swNoData <- function(beta, k, neighbors, blocks, niter=1000, random=TRUE) {
  result <- .Call( "swNoData", beta, k, neighbors, blocks, niter, random, PACKAGE = "bayesImageS")  
}

#' Fit a mixture of Gaussians to the observed data.
#' 
#' @param y A vector of observed pixel data.
#' @param niter The number of iterations of the algorithm to perform.
#' @param nburn The number of iterations to discard as burn-in.
#' @param priors A list of priors for the parameters of the model.
#' @return A matrix containing MCMC samples for the parameters of the mixture model.
#' @export
gibbsGMM <- function(y, niter=1000, nburn=500, priors=NULL) {
  result <- .Call( "gibbsGMM", y, niter, nburn, priors, PACKAGE = "bayesImageS")
}

#' Fit a univariate normal (Gaussian) distribution to the observed data.
#' 
#' @param y A vector of observed pixel data.
#' @param niter The number of iterations of the algorithm to perform.
#' @param priors A list of priors for the parameters of the model.
#' @return A list containing MCMC samples for the mean and standard deviation.
#' @examples
#' y <- rnorm(100,mean=5,sd=2)
#' res.norm <- gibbsNorm(y, priors=list(mu=0, mu.sd=1e6, sigma=1e-3, sigma.nu=1e-3))
#' summary(res.norm$mu[501:1000])
#' summary(res.norm$sigma[501:1000])
#' @export
gibbsNorm <- function(y, niter=1000, priors=NULL) {
  result <- .Call( "gibbsNorm", y, niter, priors, PACKAGE = "bayesImageS")
}

#' Calculate the sufficient statistic of the Potts model for the given labels.
#'
#' @param labels A matrix of pixel labels. 
#' @param neighbors A matrix of all neighbors in the lattice, one row per pixel.
#' @param blocks A list of pixel indices, dividing the lattice into independent blocks.
#' @param k The number of unique labels.
#' @return The sum of like neighbors.
#' @export
sufficientStat <- function(labels, neighbors, blocks, k) {
  result <- .Call( "sufficientStat", labels, neighbors, blocks, k, PACKAGE = "bayesImageS")  
}