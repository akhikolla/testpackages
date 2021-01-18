\name{EL-methods}
\docType{methods}
\alias{EL-methods}
\alias{EL}
\alias{EL,GCPM-method}
\title{Expected Loss (from Loss Distribution)}
\description{
Get the expected loss (EL) calculated from the portfolio loss distribution.
Because of the discretization and/or simulation errors, this is not equal to
the analytical EL (see \code{\link{EL.analyt}}). Please also note, that in case
of a simulative model (with Bernoulli default distribution) of the CreditRisk+
type the simulated EL tends to be smaller than the analytical one because the
conditional PD \eqn{\overline{PD}=PD\cdot (w^Tx)} has to be truncated (if \eqn{\overline{PD}>1}).
}
\usage{
EL(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{numeric value of length 1}
\seealso{\code{\link{EL.analyt}}}
\keyword{methods}
