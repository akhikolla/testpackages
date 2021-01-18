\name{GreedyMerge}
\alias{GreedyMerge}
\title{GreedyMerge}
\description{
Runs the hierarchical clustering algorithm to improve the solution obtained through GreedyICL.
}
\usage{
GreedyMerge(adj_cube, allocations, verbose = FALSE)
}
\arguments{
  \item{adj_cube}{A binary array of size \code{NxNxT} representing the dynamic network. The generic entry in position \code{[i,j,t]} is equal to \code{1} if \code{i} interacts with \code{j} during the \code{t}-th time frame, or to \code{0} otherwise. Only undirected networks with no self-edges are supported, so each slice of the array must be a symmetric matrix with null elements on the diagonal.}
  \item{allocations}{Initial allocations used by the algorithm. This should be a matrix of size \code{TxN} denoting the cluster membership of each node at each time. Values should be strictly positive integers. The value zero should be used to identify inactive nodes.}
  \item{verbose}{TRUE or FALSE, indicating whether a lenghty output should be printed out. Defaults to FALSE}
}
\value{
  \item{computing_time}{Number of seconds required to run the function.}
  \item{icl_start}{Exact log-ICL value for the initial allocations.}
  \item{icl_end}{Exact log-ICL value for the optimal allocations found.}
  \item{allocations}{Optimal allocations.}
}
\seealso{
  \code{\link{GreedyICL}}, \code{\link{CollapseLabels}}.
}
\references{
Rastelli, R. (2017) "Exact integrated completed likelihood maximisation in a stochastic block transition model for dynamic networks", \url{https://arxiv.org/abs/1710.03551}
}
\examples{
set.seed(12345)
data(reality_mining_84)
tframes <- dim(reality_mining_84)[3]
N <- dim(reality_mining_84)[1]

### Evaluate matrix containing the degrees of each node at each time
degrees <- apply(X = reality_mining_84, MARGIN = c(1,3), FUN = sum)

### create a binary matrix with ones corresponding to inactive nodes
inactive_nodes <- degrees == 0

### create an edgelist listing all of the ids (time, node_label) of inactive nodes 
n_inactive_nodes <- sum(inactive_nodes)
list_of_inactive_nodes <- matrix(NA,n_inactive_nodes,2)
index <- 1
for (t in 1:tframes) for (i in 1:N)  if (inactive_nodes[i,t])
{
  list_of_inactive_nodes[index,1] = t
  list_of_inactive_nodes[index,2] = i
  index = index + 1
}

### Find starting allocations using k-means
allocations_init <- GreedyInit(reality_mining_84, 20, list_of_inactive_nodes)

### Run the GreedyICL algorithm (this may take some time)
output_greedy <- GreedyICL(reality_mining_84, allocations_init, max_n_iter = 1)
### max_n_iter is set to 1 to speed up the demonstration: 
### please always use max_n_iter = 100 in applications

### Run the hierarchical clustering routine
output_merge <- GreedyMerge(adj_cube = reality_mining_84, allocations = output_greedy$allocations)
str(output_merge)

}



