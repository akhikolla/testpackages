\name{EdgelistFromAdjacency}
\alias{EdgelistFromAdjacency}
\title{Convert an adjacency matrix to filled edgelist.}
\usage{
  EdgelistFromAdjacency(A,
    nodelist = paste("node", 1:nrow(A), sep = ""))
}
\arguments{
  \item{A}{matrix, see 'Details' for formatting
  assumptions.}

  \item{nodelist}{character, optional list of node names.}
}
\value{
  data.frame, full list of all possible edges with weights
  for each in third column.
}
\description{
  Given the adjacency matrix for a network return a
  data.frame listing all possible edges and the weights for
  each edge.
}
\details{
  This assumes that the row of the adjacency matrix
  indicates the node the edge is coming 'from', the column
  represent the node the edge is going 'to', and the value
  in the adjacency matrix is the weight given to the edge.
}
\examples{
n <- 10
A <- matrix(rnorm(n*n), nrow=n)
A
EdgelistFromAdjacency(A)

n <- 100
A <- matrix(rnorm(n*n), nrow=n)
A
EdgelistFromAdjacency(A)

n <- 500
A <- matrix(rnorm(n*n), nrow=n)
A
\dontrun{EdgelistFromAdjacency(A)}
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/}
}
\seealso{
  \code{\link{EdgelistFromIgraph}}
}

