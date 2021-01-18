\name{coincident}
\alias{coincident}
\title{Determine coincident coordinates}
\usage{
coincident(coords1, coords2)
}
\arguments{
  \item{coords1}{An \eqn{n1 \times 2} numeric matrix of
  coordinates.}

  \item{coords2}{An \eqn{n2 \times 2} numeric matrix of
  coordinates.}
}
\value{
Returns a matrix with the indices of the coincident
locations. Specifically, an \eqn{r \times 2} matrix will be
returned, with \eqn{r} being the number of coordinates in
\code{coords1} coinciding with coordinates in
\code{coords2}. If row \code{i} of the matrix is c(2, 5),
then the \code{i}th set of coincident locations is between
the 2nd row of \code{coords1} and the 5th row of
\code{coords2}. If there are no coincident locations, then
a matrix of size \eqn{0 \times 2} is returned.
}
\description{
\code{coincident} takes the coordinate matrices
\code{coords1} and \code{coords2} and returns the indices
of the coincident coordinates of the two matrices.
}
\details{
This function calls a compiled \code{C++} program created
using the Rcpp package. This (may) result in a significant
speedup over pure R code since the search algorithm
involves loops. We assume that there are no duplicate
coordinates in \code{coords1} and \code{coords2},
respectively. In other words, each row of \code{coords1} is
unique and each row of \code{coords2} is unique. There is
at most 1 row of \code{coords1} that will match with a row
in \code{coords2}.
}
\examples{
#Generate two sets of coordinates
loc1 <- as.matrix(expand.grid(seq(0, 1, len = 25), seq(0, 1, len = 25)))
loc2 <- as.matrix(expand.grid(seq(0, 1, len = 101), seq(0, 1, len = 101)))
coincident(loc1, loc2)
}
\author{
Joshua French
}

