\name{EL.analyt-methods}
\docType{methods}
\alias{EL.analyt-methods}
\alias{EL.analyt}
\alias{EL.analyt,GCPM-method}
\title{Expected Loss (analytical)}
\description{
Get the expected loss (EL) calculated from the portfolio data. Because of the
discretization and/or simulation errors, this is not equal to the EL calculated
from the portfolio loss distribution (see \code{\link{EL}}).
}
\usage{
EL.analyt(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{numeric value of length 1}
\seealso{\code{\link{EL}}}
\keyword{methods}

