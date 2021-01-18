\name{zero}
\alias{zero}
\alias{is.zero}
\title{The zero algebraic object}
\description{
Test for a \code{freealg} object's being zero
}
\usage{
is.zero(x)
}
\arguments{
  \item{x}{Object of class \code{freealg}}
}
\details{

  Function \code{is.zero()} returns \code{TRUE} if \code{x} is indeed
  the zero free algebra object.  It is defined as
  \code{length(coeffs(x))==0} for reasons of efficiency, but
  conceptually it returns \code{x==constant(0)}.

  (Use \code{constant(0)} to create the zero object).

}
\author{Robin K. S. Hankin}
\seealso{\code{\link{constant}}}
\examples{

stopifnot(is.zero(constant(0)))

}

\keyword{symbolmath}
