\name{AdjacencyFromEdgelist}
\alias{AdjacencyFromEdgelist}
\title{Convert an edgelist to an adjacency matrix}
\usage{
  AdjacencyFromEdgelist(elist, check.full = TRUE)
}
\arguments{
  \item{elist}{data.frame, see 'Details' for formatting
  assumptions.}

  \item{check.full}{logical, if TRUE ensures that all
  possible edges are in the list exactly once; if FALSE it
  assumes this is true.}
}
\value{
  list, containing an adjacency matrix and a vector of node
  ids identifying the rows and columns. \tabular{ll}{
  \code{adjacency} \tab The adjacency matrix for the
  network. The row indicates the node the edge is coming
  'from', the column indicates the node the edge is going
  'to', and the value in the adjacency matrix is the weight
  given to the edge.\cr \code{nodelist} \tab The ids of the
  nodes in the same order as the the rows and columns of
  the adjacency matrix.\cr }
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
AdjacencyFromEdgelist(edgelist)
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

