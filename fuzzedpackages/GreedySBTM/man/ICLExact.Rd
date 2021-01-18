\name{ICLExact}
\alias{ICLExact}
\title{ICLExact}
\description{
Evaluates the exact log integrated complete likelihood for a given allocation matrix.
}
\usage{
ICLExact(adj_cube, allocations, verbose = FALSE)
}
\arguments{
  \item{adj_cube}{A binary array of size \code{NxNxT} representing the dynamic network. The generic entry in position \code{[i,j,t]} is equal to \code{1} if \code{i} interacts with \code{j} during the \code{t}-th time frame, or to \code{0} otherwise. Only undirected networks with no self-edges are supported, so each slice of the array must be a symmetric matrix with null elements on the diagonal.}
  \item{allocations}{This should be a matrix of size \code{TxN} denoting the cluster membership of each node at each time. Values should be strictly positive integers. The value zero should be used to identify inactive nodes.}
  \item{verbose}{TRUE or FALSE, indicating whether a lenghty output should be printed out. Defaults to FALSE}
}
\value{
  \item{computing_time}{Number of seconds required for the evaluation.}
  \item{prior_value}{Contribution to the exact log-ICL given by the marginal prior on the allocations.}
  \item{likelihood_value}{Contribution to the exact log-ICL given by the marginal likelihood.}
  \item{icl_value}{Exact log-ICL value for the allocation matrix given.}
}
\examples{
set.seed(12345)
data(reality_mining_84)
tframes <- dim(reality_mining_84)[3]
N <- dim(reality_mining_84)[1]
allocations <- matrix(sample(1:10,tframes*N,TRUE),tframes,N)
ICLExact(adj_cube = reality_mining_84, allocations = allocations, verbose = FALSE)
}



