\name{expSBM_ELBO}
\alias{expSBM_ELBO}
\title{expSBM_ELBO}
\description{
Evaluates the evidence lower bound for a given dynamic network.
}
\usage{
expSBM_ELBO(N, edgelist, Z, lambda, mu, nu, directed = F, trunc = T, verbose = F)
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
  \item{verbose}{\code{TRUE} or \code{FALSE} indicating whether a lengthy output should be printed out.}
}
\value{
  \item{computing_time}{Number of seconds required for the evaluation.}
  \item{elbo_value}{Value of the evidence lower bound for the given variational parameters.}
}


