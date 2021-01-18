\name{IgraphFromEdgelist}
\alias{IgraphFromEdgelist}
\title{Convert an edgelist to an igraph}
\usage{
  IgraphFromEdgelist(elist, directed = TRUE)
}
\arguments{
  \item{elist}{data.frame, see 'Details' for formatting
  assumptions.}

  \item{directed}{logical, If TRUE, the returned igraph is
  directed.}
}
\value{
  igraph, If the edgelist third column has values other
  than {0, 1} then the weights are stored in E(returned
  graph)$weight.
}
\description{
  Given the adjacency matrix for a network return a
  data.frame listing all possible edges and the weights for
  each edge.
}
\details{
  This assumes that \code{elist} is a data.frame with three
  columns. Each row is an edge in the network. The first
  column lists the node the edge is coming from, the second
  column lists the node the edge is going to, and the third
  column lists the weight of the edge.
}
\examples{
edgelist <- cbind(expand.grid(letters[1:2], letters[1:2]), runif(4))
g <- IgraphFromEdgelist(edgelist)
get.edgelist(g)
E(g)$weight
plot(g, edge.width=5*E(g)$weight, edge.curved=TRUE)
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/}
}
\seealso{
  \code{\link{EdgelistFill}}
}

