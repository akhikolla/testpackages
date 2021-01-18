\name{VaR.cont-methods}
\docType{methods}
\alias{VaR.cont-methods}
\alias{VaR.cont}
\alias{VaR.cont,GCPM-method}
\alias{VaR.cont,GCPM,numeric-method}
\title{Risk Contributions to Portfolio Value at Risk}
\description{
Get the counterparties' contributions to portfolio value at risk
(see \code{\link{GCPM-class}}). In case of a simulative model, these values are
calculated from individual losses greater or equal \code{loss.thr}
(see \code{\link{init}}). Contributions are not available if \code{loss.thr} is too high.
}
\usage{
VaR.cont(this,alpha)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
  \item{alpha}{numeric vector with entries between 0 and 1}
}
\value{numeric matrix}
\seealso{\code{\link{init}},\code{\link{loss.thr}}}
\keyword{methods}
