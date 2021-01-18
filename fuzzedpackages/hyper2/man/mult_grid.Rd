\name{mult_grid}
\alias{mult_grid}
\alias{pair_grid}
\title{
Kronecker matrix functionality
}
\description{
Peculiar version of \code{expand.grid()} for matrices
}
\usage{
mult_grid(L)
pair_grid(a,b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{List of matrices}
  \item{a,b}{Matrices}
}
\details{
Function \code{pair_grid(a,b)} returns a matrix with each column of
\code{a} \code{cbind()}-ed to each colum of \code{b}.

Function \code{mult_grid()} takes a list of matrices; it is designed for
use by \code{ggol()}.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{ggol}}}
\examples{
pair_grid(diag(2),diag(3))
mult_grid(lapply(1:4,diag))
}
