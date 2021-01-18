\name{dist1}
\alias{dist1}
\title{Calculate Euclidean distance matrix for a matrix of coordinates}
\usage{
dist1(coords)
}
\arguments{
  \item{coords}{An \eqn{nr \times nc} numeric matrix of
  coordinates.}
}
\value{
An \eqn{nr \times nr} matrix of Euclidean distances.
}
\description{
\code{dist1} takes a matrix of coordinates and returns the
Euclidean distance matrix of the coordinates.  It does this
using a compiled \code{C} program, so it is faster than the
builtin R \code{dist} function.
}
\examples{
x <- matrix(rnorm(30), ncol = 3)
dist1(x)
}
\author{
Joshua French
}
\seealso{
\link{dist}, \link{dist2}
}

