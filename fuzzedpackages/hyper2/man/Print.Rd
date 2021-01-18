\name{Print}
\alias{print.hyper2}
\alias{print}
\alias{Print}
\title{
Print methods
}
\description{
Print methods for \code{hyper2} objects
}
\usage{
\method{print}{hyper2}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{hyper2}}
  \item{\dots}{Further arguments, currently ignored}
}
\value{
Returns the \code{hyper2} object it was sent, invisibly.  Used mainly
for its side-effect 
}
\author{Robin K. S. Hankin}
\note{This is the only place that the \code{pnames} attribute is used}
\examples{
data(chess)
chess
}