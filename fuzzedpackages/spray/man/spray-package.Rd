\name{spray-package}
\alias{spray-package}
\docType{package}
\title{
Sparse arrays and multivariate polynomials
}
\description{
Functionality for sparse arrays, with emphasis on sparse arrays
interpreted as multivariate polynomials.
}
\details{

Base R has the capability of dealing with arbitrary dimensioned
numerical arrays, with the \code{array} class.
  
A sparse array is a type of array in which nonzero elements are stored
along with an index vector describing their coordinates---instead of
arrays.  This allows for efficient storage and manipulation as base
arrays often require the storing of many zero elements which consume
computational and memory resources.

One natural application for sparse arrays is multivariate polynomials
and the package vignette presents an extended discussion.

In the package, sparse arrays are represented as objects of class
\code{spray}.  They use the \code{C++} standard template library (STL)
\code{map} class, with keys being (unsigned) integer vectors, and values
floats.
}
\author{Robin K. S. Hankin}
\keyword{ package }
\examples{

# define a spray using a matrix of indices and a vector of values:
M <- matrix(sample(0:3,21,replace=TRUE),ncol=3)
a <- spray(M,sample(7))

# there are many pre-defined simple sprays:
b <- homog(3,4)

# arithmetic operators work:
a + 2*b
a - a*b^2/4
a+b

#  we can sum over particular dimensions:
asum(a+b,1)

# differentiation is supported:
deriv(a^6,2)

# extraction and replacement work as expected:

b[1,2,1]
b[1,2,1,drop=TRUE]

b[diag(3)] <- 3



}
