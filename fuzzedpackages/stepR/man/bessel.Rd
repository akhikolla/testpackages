\name{BesselPolynomial}
\alias{BesselPolynomial}

\title{Bessel Polynomials}

\description{
Recursively compute coefficients of Bessel Polynomials.

\bold{Deprecation warning:} This function is a help function for the Bessel filters in \code{\link{dfilter}} and may be removed when \code{\link{dfilter}} will be removed.
}

\usage{
BesselPolynomial(n, reverse = FALSE)
}

\arguments{
\item{n}{order}
\item{reverse}{whether to return the coefficients of a reverse Bessel Polynomial}
}

\value{
Returns the polynom's coefficients ordered increasing with the exponent, i.e. starting with the intercept, as for \code{\link{polyroot}}.
}

\seealso{\code{\link{dfilter}}, \code{\link{bessel}}, \code{\link{polyroot}}}

\examples{
# 15 x^3 + 15 x^2 + 6 x + 1
BesselPolynomial(3)
}

\keyword{math}
