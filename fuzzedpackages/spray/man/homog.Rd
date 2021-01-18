\name{homog}
\alias{homog}
\alias{product}
\alias{homog}
\alias{linear}
\alias{lone}
\alias{one}
\alias{xyz}
\title{
Various functions to create simple spray objects
}
\description{
Various functions to create simple spray objects such as single-term,
homogenous, and constant multivariate polynomials.
}
\usage{
product(power)
homog(d,power=1)
linear(x,power=1)
lone(n,d=n)
one(d)
xyz(d)
}
\arguments{
  \item{d}{An integer; generally, the dimension or arity of the resulting spray object}
  \item{power}{Integer vector of powers}
  \item{x}{Numeric vector of coefficients}
  \item{n}{In function \code{lone()}, the term to raise to power 1}
}

\value{
All functions documented here return a spray object
}
\author{
Robin K. S. Hankin
}
\note{
The functions here are related to their equivalents in the multipol
package, but are not exactly the same.

Function \code{zero()} is documented at \code{zero.Rd}, but is listed
below for convenience.

}
\seealso{\code{\link{constant}}, \code{\link{zero}}}
\examples{
product(1:3)                   #      x * y^2 * z^3
homog(3)                       #      x + y + z
homog(3,2)                     #      x^2  + xy + xz + y^2 + yz + z^2
linear(1:3)                    #      1*x + 2*y + 3*z
linear(1:3,2)                  #      1*x^2 + 2*y^2 + 3*z^2
lone(3)                        #      z
lone(2,3)                      #      y
one(3)                         #      1
zero(3)                        #      0
xyz(3)                         #      xyz
}
\keyword{symbolmath}
