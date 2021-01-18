\name{EC.cont-methods}
\docType{methods}
\alias{EC.cont-methods}
\alias{EC.cont}
\alias{EC.cont,GCPM-method}
\alias{EC.cont,GCPM,numeric-method}
\title{Risk Contributions to Economic Capital}
\description{
Calculate contributions to the economic capital on portfolio level for each
portfolio position. In case of a simulative model, the risk contributions are
calculated as contributions to expected shortfall on a lower loss level
\eqn{\tau}, such that ES(\eqn{\tau}) is as close as possible to EC(\eqn{\alpha}).
Furthermore, in case of a simulative model, loss scenarios above a predefined
threshold (\code{loss.thr}) are analyzed in order to calculate the risk
contributions. If \code{loss.thr} is too high (depending on value of \code{alpha})
the calculation will be not possible.
}
\usage{
EC.cont(this,alpha)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
  \item{alpha}{numeric vector of loss levels between 0 and 1}
}
\value{numeric matrix with number of rows equal to number of counterparties within the portfolio and number of columns equal to \code{length(alpha)}
}
\seealso{\code{\link{loss.thr}}}
\keyword{methods}
