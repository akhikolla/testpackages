\name{as.vector}
\alias{as.vector.clifford}
\alias{as.vector}
\alias{as.vector.clifford}
\title{Coerce a clifford vector to a numeric vector}
\description{
  Given a clifford object with all terms of grade 1, return the
  corresponding numeric vector
}
\usage{
\method{as.vector}{clifford}(x,mode = "any")
}
\arguments{
  \item{x}{Object of class clifford}
  \item{mode}{ignored}
}
\author{Robin K. S. Hankin}
\note{
The awkward R idiom of this function is because the terms may be stored
in any order; see the examples
}
\seealso{\code{\link{numeric_to_clifford}}}
\examples{
x <- clifford(list(6,2,9),1:3)
as.vector(x)

as.1vector(as.vector(x)) == x  # should be TRUE

}