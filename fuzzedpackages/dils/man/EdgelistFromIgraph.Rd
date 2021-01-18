\name{EdgelistFromIgraph}
\alias{EdgelistFromIgraph}
\title{Convert an igraph to filled edgelist}
\usage{
  EdgelistFromIgraph(g, useWeight = FALSE)
}
\arguments{
  \item{g}{igraph, from \link{igraph} package.}

  \item{useWeight}{logical, Should E(g)$weight be used as
  the weights for the edges?}
}
\value{
  data.frame, full list of all possible edges with weights
  for each in third column.
}
\description{
  Given an igraph object for a network return a data.frame
  listing all possible edges and the weights for each edge.
}
\details{
  This function is preferred to the igraph function
  \code{get.edgelist} because \code{get.edgelist} only
  returns rows for edges that have non-zero weight and does
  not return weights, if present.
}
\examples{
g <- erdos.renyi.game(10, 2/10)
EdgelistFromIgraph(g)

V(g)$name <- letters[1:vcount(g)]
EdgelistFromIgraph(g)

E(g)$weight <- runif(ecount(g))
EdgelistFromIgraph(g, useWeight=TRUE)
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/}
}
\seealso{
  \code{\link{EdgelistFromAdjacency}}
}

