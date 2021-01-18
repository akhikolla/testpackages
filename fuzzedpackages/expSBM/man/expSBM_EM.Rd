\name{expSBM_EM}
\alias{expSBM_EM}
\title{expSBM_EM}
\description{
Runs the variational expectation maximization algorithm for a given number of latent groups.
}
\usage{
expSBM_EM(N, edgelist, Z, lambda, mu, nu, directed = F, trunc = T, 
          tol = 0.001, n_iter_max = 100, verbose = F)
}
\arguments{
  \item{N}{Number of nodes.}
  \item{edgelist}{A matrix with 4 columns: on the first column the sender node, on the second the receiver, on the third either a one or zero to indicate whether it is an interaction or a non-interaction respectively, on the fourth the corresponding exponential length.}
  \item{Z}{A \code{NxK} matrix indicating a soft clustering of the nodes into the \code{K} latent groups. The generic entry in position \code{[i,k]} represents the posterior probability that node \code{i} belongs to group \code{k}.}
  \item{lambda}{Mixing proportions of the latent groups.}
  \item{mu}{A matrix of size \code{KxK} indicating the exponential rates for the interaction lengths, for each pair of groups. Must be a symmetric matrix if \code{directed} is false.}
  \item{nu}{A matrix of size \code{KxK} indicating the exponential rates for the non-interaction lengths, for each pair of groups. Must be a symmetric matrix if \code{directed} is false.}
  \item{directed}{\code{TRUE} or \code{FALSE} indicating whether interactions have an orientation or not.}
  \item{trunc}{\code{TRUE} or \code{FALSE} indicating whether the first and last interactions or non-interactions for every edge are assumed to be truncated or not.}
  \item{tol}{Stop the maximization if the relative increase in the objective function is not larger than this value.}
  \item{n_iter_max}{Stop the maximization if the number of iterations is larger than this value. This parameter can be set to zero or one for debug purposes.}
  \item{verbose}{\code{TRUE} or \code{FALSE} indicating whether a lengthy output should be printed out.}
}
\value{
  \item{computing_time}{Number of seconds required for the evaluation.}
  \item{elbo_values}{Stored values of the objective function at each iteration.}
  \item{Z_star}{Optimal soft clustering of the nodes into the groups.}
  \item{lambda_star}{Optimal mixing proportions.}
  \item{mu_star}{Optimal group-specific parameters for the exponential rates of the interaction lengths.}
  \item{nu_star}{Optimal group-specific parameters for the exponential rates of the non-interaction lengths.}
}
\examples{
set.seed(1)
data(high_school)
K <- 4
lambda_init <- rep(1/K,K)
Z_init <- expSBM_init(high_school$edgelist, K, "random")
mu_init <- nu_init <- matrix(1,K,K)
expSBM_EM(N = 327, high_school$edgelist, Z_init, lambda_init, mu_init, nu_init)
}
