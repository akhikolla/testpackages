\name{ES.cont-methods}
\docType{methods}
\alias{ES.cont-methods}
\alias{ES.cont}
\alias{ES.cont,GCPM-method}
\alias{ES.cont,GCPM,numeric-method}
\title{Risk Contributions to Expected Shortfall}
\description{
Calculate contributions to the expected shortfall on portfolio level for each
portfolio position. In case of a simulative model, loss scenarios above a
predefined threshold (\code{loss.thr}) are analyzed in order to calculate the
risk contributions. If \code{loss.thr} is too high, calculation may be not
possible (depending on value of \code{alpha}). 
}
\usage{
ES.cont(this,alpha)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
  \item{alpha}{numeric vector of loss levels between 0 and 1}
}
\value{numeric matrix with number of rows equal to number of counterparties within the portfolio and number of columns equal to \code{length(alpha)}
}
\seealso{\code{\link{loss.thr}}}
\keyword{methods}
