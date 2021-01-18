\name{GreedyICL}
\alias{GreedyICL}
\title{GreedyICL}
\description{
Runs the GreedyICL algorithm on the provided network to maximise the exact integrated complete likelihood.
}
\usage{
GreedyICL(adj_cube, allocations, max_n_iter = 100, verbose = FALSE)
}
\arguments{
  \item{adj_cube}{A binary array of size \code{NxNxT} representing the dynamic network. The generic entry in position \code{[i,j,t]} is equal to \code{1} if \code{i} interacts with \code{j} during the \code{t}-th time frame, or to \code{0} otherwise. Only undirected networks with no self-edges are supported, so each slice of the array must be a symmetric matrix with null elements on the diagonal.}
  \item{allocations}{Initial allocations used by the greedy algorithm. This should be a matrix of size \code{TxN} denoting the cluster membership of each node at each time. Values should be strictly positive integers. The value zero should be used to identify inactive nodes.}
  \item{max_n_iter}{A positive integer denoting the maximum number of iterations for the GreedyICL maximisation. Defaults to \code{100} (i.e. each node is updated up to \code{100} times).}
  \item{verbose}{TRUE or FALSE, indicating whether a lenghty output should be printed out. Defaults to FALSE}
}
\value{
  \item{computing_time}{Number of seconds required to perform the optimisation.}
  \item{icl_start}{Exact log-ICL value for the initial allocations.}
  \item{icl_trace}{Exact log-ICL values for each of the solutions visited during the optimisation.}
  \item{icl_end}{Exact log-ICL value for the optimal allocations found.}
  \item{allocations}{Optimal allocations.}
}
\seealso{
  \code{\link{GreedyInit}}, \code{\link{GreedyMerge}}, \code{\link{CollapseLabels}}.
}
\references{
Rastelli, R. (2017) "Exact integrated completed likelihood maximisation in a stochastic block transition model for dynamic networks", \url{https://arxiv.org/abs/1710.03551}
}
\examples{
### A complete example is provided in the documentation of GreedyMerge.
}



