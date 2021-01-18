\name{dist2}
\alias{dist2}
\title{Calculate Euclidean distance matrix between coordinates of two matrices}
\usage{
dist2(coords1, coords2)
}
\arguments{
  \item{coords1}{An \eqn{nr1 \times nc1} numeric matrix of
  coordinates.}

  \item{coords2}{An \eqn{nr2 \times nc2} numeric matrix of
  coordinates.}
}
\value{
An \eqn{nr1 \times nr2} matrix of Euclidean distances.
}
\description{
\code{dist2} takes the matrices of coordinates
\code{coords1} and \code{coords2} and returns the
inter-Euclidean distances between coordinates.
}
\examples{
x1 <- matrix(rnorm(30), ncol = 3)
x2 <- matrix(rnorm(60), ncol = 3)
dist2(x1, x2)
}
\author{
Joshua French
}
\seealso{
dist, dist1
}

