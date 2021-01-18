\name{GenerateDilsNetwork}
\alias{GenerateDilsNetwork}
\title{Combine multiple networks into a single weighted network.}
\usage{
  GenerateDilsNetwork(x, subsample = 10000,
    n.subsamples = 1000, ignore.cols, use.cols,
    progress.bar = FALSE)
}
\arguments{
  \item{x}{data.frame, data over which to run PCA}

  \item{subsample}{numeric or logical, If an integer, size
  of each subsample.  If FALSE, runs PCA on entire data
  set.}

  \item{n.subsamples}{numeric, number of subsamples.}

  \item{ignore.cols}{numeric, indices of columns not to
  include}

  \item{use.cols}{numeric, indices of columns to use}

  \item{progress.bar}{logical, if TRUE then progress in
  running subsamples will be shown.}
}
\value{
  A list \item{dils}{vector, named vector of component
  weights for first dimension of principal component
  analysis (see example for comparison to
  \code{\link{prcomp}}).} \item{dils.edgelist}{Unused
  columns of \code{x} bound with the DILS scores on the
  right. Forms an edgelist if there were two unused columns
  and they containted the ids for the dyads.}
  \item{coefficients}{named vector, weights that genereate
  \code{dils} by taking dot-product with network data.}
  \item{weights}{named vector, raw.weights scaled by
  standard deviations of network edges, then scaled to sum
  to 1.}
}
\description{
  Use ScalablePCA to recover optimal weights for each
  network, then calculate the weighted average across
  networks for each edge.
}
\examples{
data(iris)        # provides example data
GenerateDilsNetwork(iris, subsample=10, use.cols=1:4)
GenerateDilsNetwork(iris, subsample=10, ignore.cols=5)
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{https://github.com/shaptonstahl/}
}
\seealso{
  \code{\link{prcomp}}
}

