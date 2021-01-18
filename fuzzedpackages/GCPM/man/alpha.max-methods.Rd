\name{alpha.max-methods}
\docType{methods}
\alias{alpha.max-methods}
\alias{alpha.max,GCPM}
\alias{alpha.max}
\title{Maximum CDF Level}
\description{
Get the maximum value of the model's CDF. For simulative models,
the value should be equal to 1. For an analytical model,
the value depends on the value specified during initiation of the model (see \code{\link{init}}).
}
\usage{
  alpha.max(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{numeric of length 1}
\keyword{methods}
\seealso{\code{\link{init}}}

