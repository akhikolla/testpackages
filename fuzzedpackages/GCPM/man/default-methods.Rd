\name{default-methods}
\docType{methods}
\alias{default-methods}
\alias{default}
\alias{default,GCPM-method}
\title{Default Distribution}
\description{
Get the default distribution of each portfolio position. Using \dQuote{Poisson}
as default distribution one can simulate the standard CR+ model or group
smaller counterparties into a pool and simulate their defaults.
}
\usage{
  default(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{character of length equal to number of portfolio positions}
\seealso{\code{\link{portfolio.pois}}}
\keyword{methods}
