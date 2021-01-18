\name{tt.events}
\alias{tt.events}
\docType{data}
\title{
    Data.frame of hosts events information by time.
    Fake data.
}
\description{
  Fake simulated data of hosts events over the time.
}
\usage{data("transmissiontree")}
\format{
  A data frame with 63 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{a character vector}
    \item{\code{time}}{a character vector}
    \item{\code{status}}{a character vector}
    \item{\code{infectedby}}{a character vector}
    \item{\code{probabilities}}{a character vector}
    \item{\code{X}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
  }
}
\examples{
data(transmissiontree)
print(tt.events)
}
\keyword{datasets}
