\name{LHR-methods}
\docType{methods}
\alias{LHR-methods}
\alias{LHR}
\alias{LHR,GCPM-method}
\title{Likelihood Ratio}
\description{
Get the likelihood ratio for each scenario defined in \code{random.numbers} 
(see \code{\link{init}})}
\usage{
LHR(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{numeric vector of length equal to \code{nrow(random.numbers)}}
\keyword{methods}
