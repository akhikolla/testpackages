\name{constant}
\alias{constant}
\alias{constant}
\alias{is.constant}
\alias{constant<-}
\alias{constant.freealg}
\alias{constant<-.freealg}
\alias{constant.numeric}
\title{The constant term}
\description{Get and set the constant term of a \code{freealg} object}
\usage{
\method{constant}{freealg}(x)
\method{constant}{numeric}(x)
\method{constant}{freealg}(x) <- value
is.constant(x)
}
\arguments{
  \item{x}{Object of class \code{freealg}}
  \item{value}{Scalar value for the constant}
}
\details{

  The constant term in a free algebra object is the coefficient of the
  empty term.  In a \code{freealg} object, the map including \code{\{\}
  -> v} implies that \code{v} is the constant.

  If \code{x} is a \code{freealg} object, \code{constant(x)} returns
  the value of the constant in the multivariate polynomial; if \code{x}
  is numeric, it returns a constant \code{freealg} object with value
  \code{x}.

  Function \code{is.constant()} returns \code{TRUE} if its argument has
  no variables and \code{FALSE} otherwise.
  
}
\author{Robin K. S. Hankin}
\examples{
p <- as.freealg("1+X+Y+xy")

constant(p)
constant(p^5)

constant(p) <- 1000
p


}
