\name{MeasureNetworkInformation}
\alias{MeasureNetworkInformation}
\title{Measure how much a network informs a particular network measure}
\usage{
  MeasureNetworkInformation(g, FUN = betweenness,
    remove.share = 0.2, sample.size = 100,
    progress.bar = FALSE)
}
\arguments{
  \item{g}{igraph, graph to measure}

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
  numeric, mean precision of the measure \code{FUN} across
  the network
}
\description{
  Given an \code{igraph} network, repeatedly perturb the
  graph and take some measure of the network to see how
  much the measure varies, then return a measure that
  increases as the precision of the function values
  increases.
}
\details{
  This function can vary tremendously based on the network
  measure being considered and the other parameters.  It is
  only recommended that this be used for comparing the
  informativeness of two networks on the same set of nodes,
  keeping all the parameters the same.

  Here information is measured as 1 / mean across and
  perturbed graphs nodes of the relative error of a network
  node measure.

  Specifically, \code{FUN} is applied to the graph \code{g}
  to generate reference values.  Some \code{sample.size}
  copies of the igraph are generated.  For each,
  \code{round(remove.share * n.edges)} randomly selected
  edges are dropped to generate a perturbed graph.  For
  each perturbed graph \code{FUN} is applied, generating a
  value for each node in the network.  For each node the
  relative error

  \deqn{|\frac{measure of perturbed g - measure of
  g}{measure of g}|}{abs( (measure of perturbed g - measure
  of g) / measure of g )}

  is calculated, then the mean of these across nodes and
  perturbed graphs is calculated, generating a mean
  relative error for the network. This value is
  reciprocated to get a measure of precision.

  This measure appears to be very sensitive to the choice
  of \code{FUN}.
}
\examples{
g.rand <- random.graph.game(100, 5/100)
m.rand <- MeasureNetworkInformation(g.rand)
m.rand

pf <- matrix( c(.8, .2, .3, .7), nr=2)
g.pref <- preference.game(100, 2, pref.matrix=pf)
m.pref <- MeasureNetworkInformation(g.pref)
m.pref

m.pref / m.rand  # Relative informativeness of this preference graph
                 # to this random graph with respect to betweenness
\dontrun{
prob.of.link <- c(1:50)/100
mnis <- sapply(prob.of.link, function(p)
  MeasureNetworkInformation(random.graph.game(100, p)))
plot(prob.of.link, mnis,
     type="l",
     main="Network Information of random graphs",
     xlab="probability of link formation",
     ylab="information")
mtext("with respect to betweenness measure", line=0.5)}
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/dils}
}

