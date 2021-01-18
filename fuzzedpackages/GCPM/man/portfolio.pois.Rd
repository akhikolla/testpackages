\name{portfolio.pois}
\alias{portfolio.pois}
\docType{data}
\title{
Example Portfolio Data with Poisson Default Mode}
\description{
The dataset contains an example portfolio in the structure needed by the \code{\link{analyze}} function.}
\usage{data("portfolio.pois")}
\format{
  A data frame with 3000 counterparties and the following variables.
  \describe{
    \item{\code{Number}}{Counterparty ID (numeric)}
    \item{\code{Name}}{Counterparty name (character)}
    \item{\code{Business}}{Business line (character)}
    \item{\code{Country}}{Country (character)}
    \item{\code{EAD}}{Exposure at default (numeric)}
    \item{\code{LGD}}{Loss given dafault (numeric)}
    \item{\code{PD}}{Probability of default (numeric)}
    \item{\code{Default}}{Default mode (\sQuote{Poisson} or \sQuote{Benroulli})}
    \item{\code{A}}{sector weights for sector A}
    \item{\code{B}}{sector weights for sector B}
    \item{\code{C}}{sector weights for sector C}
  }
}
\keyword{datasets}
