\name{ooom}
\alias{ooom}
\title{
One-over-one-minus for spray objects
}
\description{
One-over-one-minus for spray objects; the nearest to \sQuote{division}
that we can get.
}
\usage{
ooom(S, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{object of class spray}
  \item{n}{Order of the approximation}
}
\details{
Returns the Taylor expansion to order \eqn{n}.
}
\value{
Returns a spray object of the same arity as \code{S}.
}
\author{
Robin K. S. Hankin
}
\note{
Uses Horner's method for efficiency
}
\examples{

a <- homog(4,2)
jj <-  (1-a)*ooom(a,3)


constant(jj)    # should be 1
rowSums(index(jj))   # a single 0 and lots of 8s.


}

\keyword{symbolmath}
