\name{RelativeNetworkInformation}
\alias{RelativeNetworkInformation}
\title{Compare how much two networks inform a particular network measure}
\usage{
  RelativeNetworkInformation(g1, g2, FUN = betweenness,
    remove.share = 0.2, sample.size = 100,
    progress.bar = FALSE)
}
\arguments{
  \item{g1}{igraph, graph to measure}

  \item{g2}{igraph, graph to measure}

  \item{FUN}{function, a function that takes an igraph and
  returns a value for each node in the network.}

  \item{remove.share}{numeric, fraction of the edges that
  are removed randomly when perturbing the network.}

  \item{sample.size}{numeric, number of perturbed graphs to
  generate}

  \item{progress.bar}{logical, if TRUE then a progress bar
  is shown.}
}
\value{
  list, containing the following \tabular{lll}{
  \code{g1.over.g2} \tab numeric \tab informativeness of
  the first network over the second \cr \code{winner} \tab
  character \tab either g1 or g2 \cr \code{g1.measure} \tab
  numeric \tab \code{MeasureNetworkInformation(g1, ...)}
  \cr \code{g2.measure} \tab numeric \tab
  \code{MeasureNetworkInformation(g2, ...)} \cr }
}
\description{
  Given two \code{igraph} networks, use
  \code{\link{MeasureNetworkInformation}} to gauge the
  informativeness of each network and then return the
  ratio. If greater than 1, then the first network
  specified is more informative.
}
\details{
  This measure appears to be very sensitive to the choice
  of \code{FUN}. See
  \code{\link{MeasureNetworkInformation}} for details.
}
\examples{
g.rand <- random.graph.game(100, 5/100)

pf <- matrix( c(.8, .2, .3, .7), nr=2)
g.pref <- preference.game(100, 2, pref.matrix=pf)

RelativeNetworkInformation(g.rand, g.pref)
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/dils}
}

