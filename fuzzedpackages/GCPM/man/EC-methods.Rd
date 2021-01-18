\name{EC-methods}
\docType{methods}
\alias{EC-methods}
\alias{EC}
\alias{EC,GCPM,missing-method}
\alias{EC,GCPM,numeric-method}
\title{Economic Capital}
\description{Get the value of economic capital for the portfolio on level(s)
\code{alpha}}
\usage{
EC(this,alpha)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
  \item{alpha}{numeric vector of loss levels between 0 and 1}
}
\value{numeric vector of length equal to \code{length(alpha)}.
}
\keyword{methods}

