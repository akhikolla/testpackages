\name{regmed_example}
\alias{regmed_example}
\docType{data}
\title{
Example dataset for regmed package
}
\description{
Example methylation Dataset from public data source with 85 observations, and narrowed down
to 6 variables for the example.
}
\usage{data("regmed_example")}
\format{
  A data frame with 85 observations on the following 6 variables.
  \describe{
    \item{\code{x}}{a numeric vector for exposure}
    \item{\code{y}}{a numeric vector for stress response}
    \item{\code{med.cg01644731}}{a numeric vector for methylation as a mediator}
    \item{\code{med.cg06890779}}{a numeric vector for methylation as a mediator}
    \item{\code{med.cg00000029}}{a numeric vector for methylation as a mediator}
    \item{\code{med.cg00000108}}{a numeric vector, for methylation as a mediator}
  }
}
\source{
  E-GEOD-77445 - Genome-wide DNA methylation levels.
}
\keyword{datasets}
