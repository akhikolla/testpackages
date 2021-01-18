\name{idiosyncr-methods}
\docType{methods}
\alias{idiosyncr-methods}
\alias{idiosyncr}
\alias{idiosyncr,GCPM-method}
\title{Idiosyncratic Risk Weights}
\description{
Get the idiosyncratic risk weights (i.e. risk weights which are not assigned
to any sector). Currently only available if \code{model.type} == "CRP".}
\usage{
idiosyncr(this)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
}
\value{numeric vector of length equal to number of counterparties}
\keyword{methods}
