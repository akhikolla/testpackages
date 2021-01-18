\name{SimulatedData}
\alias{SimulatedData}
\docType{data}
\title{
Simulated data
}
\description{
A simulated multivariate longitudinal data for demonstration.
}
\usage{data("SimulatedData")}
\format{
  A \code{list} consists of \code{Y} the observations 100 subjects in 3 attributes along 10 time points, 
  \code{X} the design matrix with 6 covariate including the intercept, \code{TimePointsAvailable} 
  the avilable time points for each subject. 
  \describe{
    \item{\code{Y}}{The response variables.}
    \item{\code{X}}{The design matrix.}
    \item{\code{TimePointsAvailable}}{The available time points for each subject.}
  }
}

\examples{
library(MLModelSelection)
data(SimulatedData)
SimulatedData = data(SimulatedData)
}
\keyword{datasets}
