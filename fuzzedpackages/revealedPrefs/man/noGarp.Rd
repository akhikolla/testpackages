\name{noGarp}
\alias{noGarp}
\alias{noWarp}
\alias{noSarp}
\alias{okSarp}
\alias{noAxiom}
\docType{data}
\title{revealedPrefs example datasets}

\description{Toy examples to test for consistency with rationality axioms.}

\usage{
data(noGarp)
data(noWarp)
data(noSarp)
data(okSarp)
data(noAxiom)
}

\format{
Each dataset is a list of two matrices, named x and p, that contain quantities and prices, respectively. Each row is an observation, each column is a type of good.
}

\details{
Data in \code{noAxiom} violate WARP, GARP and SARP.

Data in \code{noWarp} violate WARP and SARP, but not GARP.

Data in \code{noGarp} violate GARP and SARP, but not WARP.

Data in \code{noSarp} violate SARP, but not WARP or GARP.

Data in \code{okSarp} are consistent with WARP, SARP and GARP.
}
\examples{
data(noWarp)
checkWarp(noWarp$x, noWarp$p)
}
\author{
Julien Boelaert \email{jubo.stats@gmail.com}
}

\keyword{datasets}
