\name{VaR-methods}
\docType{methods}
\alias{VaR-methods}
\alias{VaR}
\alias{VaR,GCPM-method}
\alias{VaR,GCPM,numeric-method}
\alias{VaR,GCPM,missing-method}
\title{Portfolio Value at Risk}
\description{
Calculate the portfolio value at risk on level(s) \code{alpha}.
}
\usage{
VaR(this,alpha)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
  \item{alpha}{numeric vector with entries between 0 and 1}
}
\value{numeric value of length equal to length of \code{alpha}}
\keyword{methods}
