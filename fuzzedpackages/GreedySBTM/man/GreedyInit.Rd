\name{GreedyInit}
\alias{GreedyInit}
\title{GreedyInit}
\description{
Initialises the allocation variables using k-means.
}
\usage{
GreedyInit(adj_cube, Kup, list_of_inactive_nodes = NULL)
}
\arguments{
  \item{adj_cube}{A binary array of size \code{NxNxT} representing the dynamic network. The generic entry in position \code{[i,j,t]} is equal to \code{1} if \code{i} interacts with \code{j} during the \code{t}-th time frame, or to \code{0} otherwise. Only undirected networks with no self-edges are supported, so each slice of the array must be a symmetric matrix with null elements on the diagonal.}
  \item{Kup}{Number of groups: must be a positive integer between 1 and N.}
  \item{list_of_inactive_nodes}{A matrix whose rows identify all of the inactive nodes in the network. The first element of each row denotes the time frame \code{t}, while the second element the node label \code{i}: the inactive nodes are identified by all such pairs \code{(t,i)}.}
}
\value{
The function returns a \code{TxN} allocation matrix, with zeros corresponding to inactive nodes.
}
\seealso{
  \code{\link{GreedyICL}}, \code{\link{GreedyMerge}}, \code{\link{CollapseLabels}}.
}
\examples{
### A complete example is provided in the documentation of GreedyMerge.
}



