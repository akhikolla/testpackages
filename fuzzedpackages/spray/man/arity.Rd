\name{arity}
\alias{arity}
\title{
The arity of a spray object
}
\description{
The arity of a spray object: the number of indices needed to retrieve an
entry, or the number of columns in the index matrix.
}
\usage{
arity(S)
}
\arguments{
  \item{S}{a spray object}
}
\value{
Returns an integer
}
\author{
Robin K. S. Hankin
}
\examples{
arity(spray(diag(1:6)))

}
\keyword{symbolmath}

