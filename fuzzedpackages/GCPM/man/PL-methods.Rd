\name{PL-methods}
\docType{methods}
\alias{PL-methods}
\alias{PL}
\alias{PL,GCPM-method}
\title{Counterparty Potential Loss}
\description{
Get the value of \code{PL}, the potential losses of counterparties 
(see \code{\link{GCPM-class}}). Please note, that the potential losses are 
discretized according to \code{loss.unit} (see \code{\link{init}}).
}
\usage{
PL(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{numeric value of length equal to the number of counterparties}
\seealso{\code{\link{portfolio.pois}},\code{\link{init}}}
\keyword{methods}