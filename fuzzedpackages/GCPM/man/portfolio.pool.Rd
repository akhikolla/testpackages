\name{portfolio.pool}
\alias{portfolio.pool}
\docType{data}
\title{
Pooled Portfolio}
\description{
In order to speed up calculations, counterparties of \code{\link{portfolio.pois}} with EAD*LGD < 200,000 are grouped together (pooled).
}
\usage{data("portfolio.pool")}
\format{
   A data frame with 1400 counterparties and 3 pools (each per sector) and the following variables.
  \describe{
    \item{\code{Number}}{Counterparty ID (numeric)}
    \item{\code{Name}}{Counterparty name (character)}
    \item{\code{Business}}{Business line (character)}
    \item{\code{Country}}{Country (character)}
    \item{\code{EAD}}{Exposure at default (numeric); pool: average EAD per counterparty}
    \item{\code{LGD}}{Loss given dafault (numeric); pool: EAD-weighted average LGD per counterparty}
    \item{\code{PD}}{Probability of default (numeric); pool: expectation of number of defaults}
    \item{\code{Default}}{Default mode (\sQuote{Poisson} for pools or \sQuote{Benroulli})}
    \item{\code{A}}{sector weights for sector A}
    \item{\code{B}}{sector weights for sector B}
    \item{\code{C}}{sector weights for sector C}
  }
}
\keyword{datasets}
