\name{slpm_gen}
\alias{slpm_gen}
\title{slpm_gen}
\description{
Generates the adjacency matrix \code{adj} of a SparseLPM by sampling both the data and model parameters from the posterior distribution.
}
\usage{
slpm_gen(M, N, K, hyper_pars = NULL)
}
\arguments{
  \item{M}{Number of rows of \code{adj}.}
  \item{N}{Number of cols of \code{adj}.}
  \item{K}{Number of latent dimensions of the SparseLPM.}
  \item{hyper_pars}{A list defining the hyperparameters of the model. If left as \code{NULL} all the hyperparameters are initialised to \code{1}. Otherwise, the list should contain three vectors of \code{K} positive values denoted \code{delta}, \code{a_gamma} and \code{b_gamma}, respectively.}
}
\value{
A list with components:
  \item{adj}{Generated adjacency matrix.}
  \item{U}{Generated latent positions for senders.}
  \item{V}{Generated latent positions for receivers.}
  \item{lambda}{Latent variables attached to each of the edges, selecting which dimension determines the edge probability.}
  \item{gamma}{Vector of the Gaussian precisions associated to the latent dimensions.}
}
\examples{
set.seed(12345)
network <- slpm_gen(M = 10, N = 8, K = 2)
}


