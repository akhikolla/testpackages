\name{print}
\alias{print.clifford}
\alias{as.character}
\alias{as.character.clifford}
\alias{catterm}
\title{Print clifford objects}
\description{Print methods for Clifford algebra}
\usage{
\method{print}{clifford}(x,...)
\method{as.character}{clifford}(x,...)
catterm(a)
}
\arguments{
  \item{x}{Object of class \code{clifford} in the print method}
  \item{...}{Further arguments, currently ignored}
  \item{a}{Integer vector representing a term}
}
\author{Robin K. S. Hankin}
\note{

The print method does not change the internal representation of a
\code{clifford} object, which is a two-element list, the first of which
is a list of integer vectors representing terms, and the second is a
numeric vector of coefficients.

The print method has special dispensation for length-zero clifford
objects.  It is sensitive to the value of \code{options("separate")}
which, if \code{TRUE} prints the basis vectors separately and otherwise
prints in a compact form.  The indices of the basis vectors are
separated with \code{option("basissep")} which is usually \code{NULL}
but if \eqn{n>9}, then setting \code{options("basissep" = ",")} might
look good.

Function \code{as.character.clifford()} is also sensitive to these
options.

Function \code{catterm()} is a low-level helper function.


}
\seealso{\code{\link{clifford}}}  
\examples{

rcliff()   # fine
rcliff(d=15) # incomprehensible


options("separate" = TRUE)
rcliff(d=15) # incomprehensible


options("separate" = FALSE)


}
