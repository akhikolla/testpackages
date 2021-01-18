\name{fillup}
\alias{fillup}
\alias{indep}
\title{Fillup function}
\description{
  Function \code{fillup()} concatenates a vector with a \sQuote{fillup}
  value to ensure a unit sum; if given a matrix, attaches a column so
  the rowsums are 1.

  Function \code{indep()} is the inverse: it removes the final element
  of a vector, leaving only an independent set.
  
}
\usage{
fillup(x)
indep(x)
}
\arguments{
  \item{x}{Numeric vector}
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{equalp}}}
\examples{
fillup(c(1/2,1/3))

indep(c(1/2,1/3,1/6))

}