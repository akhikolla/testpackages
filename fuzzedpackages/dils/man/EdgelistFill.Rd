\name{EdgelistFill}
\alias{EdgelistFill}
\title{Ensure an edgelist has all dyads and a column of weights.}
\usage{
  EdgelistFill(elist, fillBlanksWith = 0, nodelist)
}
\arguments{
  \item{elist}{data.frame or matrix, see 'Details' for
  formatting assumptions.}

  \item{fillBlanksWith}{numeric, default weight for edges
  not already listed in elist.}

  \item{nodelist}{character, optional list of node names.}
}
\value{
  data.frame, full list of all possible edges with weights
  for each in third column.
}
\description{
  Given a matrix or data.frame edgelist, fill in all
  possible edges not already listed with a weight of 0 or
  the value of \code{fillBlanksWith}.
}
\details{
  The \code{elist} can be either a data.frame or a matrix
  with either 2 or 3 columns. Each row is an edge. The
  first column lists the node the edge is 'from' and the
  second column lists the node the edge is 'to'. If there
  is a third column, it lists the weight of the edge.
}
\examples{
g <- erdos.renyi.game(10, 2/10)
EdgelistFill(get.edgelist(g))
EdgelistFill(get.edgelist(g), nodelist=1:10)

E(g)$weight <- runif(ecount(g))
el <- cbind(get.edgelist(g), E(g)$weight)
EdgelistFill(el)
EdgelistFill(el, nodelist=1:10)
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/}
}

