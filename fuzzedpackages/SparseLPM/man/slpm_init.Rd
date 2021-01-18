\name{slpm_init}
\alias{slpm_init}
\title{slpm_init}
\description{
Initialises the variational parameters of a SparseLPM.
}
\usage{
slpm_init(X, K, method = "random", threshold = 0.1, stdev = NULL)
}
\arguments{
  \item{X}{Rectangular adjacency matrix with non-negative entries.}
  \item{K}{Number of latent dimensions of the SparseLPM.}
  \item{method}{The variational parameters are initialised at random. However, if \code{method="distance"}, a distance-based method is used as described in Rastelli ... (2018).}
  \item{threshold}{A small number added to each of the entries of \code{X} to avoid numerical instability.}
  \item{stdev}{Standard deviation of the Gaussian used to generate the random latent positions.}
}
\value{Returns a list of variational parameters that can be used as input for \link{slpm_nga} or \link{slpm_elbo}.
}
\examples{
set.seed(12345)
M <- N <- 10
K <- 2
network <- slpm_gen(M = M, N = N, K = K)
var_pars_init <- slpm_init(X = network$adj, K = K)
}


