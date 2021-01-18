\name{SD.cont-methods}
\docType{methods}
\alias{SD.cont-methods}
\alias{SD.cont}
\alias{SD.cont,GCPM-method}
\title{Risk Contributions to Portfolio Standard Deviation}
\description{
Get the counterparties' contributions to portfolio standard deviation
(see \code{\link{GCPM-class}}). These values are only available in case of an
analytical model. 
}
\usage{
SD.cont(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{numeric value of length equal to number of counterparties}
\keyword{methods}