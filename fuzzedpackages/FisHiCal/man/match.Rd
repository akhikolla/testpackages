\name{match}
\alias{match}
\docType{data}
\title{
in-silico FISH distances and matching Hi-C frequencies with 5\% noise
}
\description{
in-silico distances were computed from the random configuration  \code{\link{conf}}. A power law model was then used to generate matching in-silico Hi-C frequencies. Random Noise was further added to long range frequencies to mimic a typical situation for Hi-C data. This data structure is used to illustrate the usage of \code{\link{prepareCalib}} and \code{\link{calibrate}}.
}
\usage{data(match)}
\format{
  A data frame with 4950 observations on the following 2 variables.
  \describe{
  \item{\code{distances}}{a numeric vector giving the in-silico FISH distances}.
    \item{\code{frequencies}}{a numeric vector giving the in-silico Hi-C contact frequencies}.
   }
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\examples{
data(match)
plot(match$frequencies ~ match$distances, xlab = "distances", ylab = "frequencies")
}
\keyword{datasets}
\keyword{datasets}
