\name{ES-methods}
\docType{methods}
\alias{ES-methods}
\alias{ES}
\alias{ES,GCPM,missing-method}
\alias{ES,GCPM,numeric-method}
\title{Expected Shortfall}
\description{Get the value of the expected shortfall for the portfolio on level(s)
\code{alpha}}
\usage{
ES(this,alpha)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
  \item{alpha}{numeric vector of loss levels between 0 and 1}
}
\value{numeric vector of length equal to \code{length(alpha)}.
}
\keyword{methods}
