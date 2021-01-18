\name{print.spray}
\alias{print.spray}
\alias{sprayvars}
\alias{print_spray_matrixform}
\alias{print_spray_polyform}
\title{
Print methods for spray objects
}
\description{
Print methods for spray objects with options for printing in matrix form
or multivariate polynomial form
}
\usage{
\method{print}{spray}(x, ...)
print_spray_matrixform(S)
print_spray_polyform(S)
}
\arguments{
  \item{x,S}{spray object}
  \item{...}{Further arguments (currently ignored)}
}
\details{

 The print method, \code{print.spray()}, dispatches to helper functions
  \code{print_spray_matrixform()} and \code{print_spray_polyform()}
  depending on the value of option \code{polyform}; see the examples
  section.

  Option \code{sprayvars} is a character vector with entries
  corresponding to the variable names for printing.

  Note that printing a spray object (in either matrix form or polynomial
  form) generally takes much longer than calculating it.
}
\value{
Returns its argument invisibly.
}
\author{
Robin K. S. Hankin
}
\note{
  There are a couple of hard-wired symbols for multiplication and
  equality which are defined near the top of the helper functions.
  
}
\examples{

(a <- spray(diag(3)))


options(polyform = FALSE)
a^3

options(polyform = TRUE)
a^3


options(sprayvars=letters)
a <- diag(26)
spray(a)

## Following example from mpoly:
a[1 + cbind(0:25, 1:26) \%\% 26] <- 2
spray(a)



}
\keyword{symbolmath}
