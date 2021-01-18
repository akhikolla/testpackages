\name{PD-methods}
\docType{methods}
\alias{PD-methods}
\alias{PD}
\alias{PD,GCPM-method}
\title{Counterparty Probability of Default}
\description{
Get the value of \code{PD}, the counterparties default probabilities within the 
portfolio (see \code{\link{analyze}}. Please note, that these PDs are adjusted 
because of discretization in order to preserve the expected loss.)
}
\usage{
PD(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{numeric value of length equal to the number of counterparties}
\seealso{\code{\link{portfolio.pois}}}
\keyword{methods}