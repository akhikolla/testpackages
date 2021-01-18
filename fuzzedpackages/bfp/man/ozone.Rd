\name{ozone}
\alias{ozone}
\docType{data}
\title{Ozone data from Breiman and Friedman, 1985}
\description{
This is the Ozone data discussed in Breiman and Friedman (JASA, 1985, p. 580).
These data are for 330 days in 1976.  All measurements are in the area of
Upland, CA, east of Los Angeles.
}
\format{
  A data frame with 366 observations on the following 13 variables.
  \describe{
    \item{\code{month}}{month of the year}
    \item{\code{day}}{day of the month}
    \item{\code{weekday}}{day of the week: a factor with levels \code{Monday},
      \code{Tueday}, \code{Wednesday}, \code{Thursday}, \code{Friday},
      \code{Saturday}, \code{Sunday}} 
    \item{\code{hourAverageMax}}{maximum 1-hour average ozone level [ppm]}
    \item{\code{pressure500Height}}{500 millibar pressure height [meters]}
    \item{\code{windSpeed}}{wind speed [mph]}
    \item{\code{humidity}}{relative humidity [\%]}
    \item{\code{tempSandburg}}{temperature at Sandberg, CA [degrees F]}
    \item{\code{tempElMonte}}{temperature at El Monte, CA [degrees F]}
    \item{\code{inversionBaseHeight}}{inversion base height [feet]}
    \item{\code{pressureGradientDaggett}}{pressure gradient from LAX to
      Daggett, CA [mm Hg]}
    \item{\code{inversionBaseTemp}}{inversion base temperature [degrees F]}
    \item{\code{visibility}}{visibility [miles]}
  }
}


\usage{data(ozone)}
\source{Breiman, L and  Friedman, J.  (1985),
 \dQuote{Estimating Optimal Transformations for Multiple Regression and Correlation},
 \emph{Journal of the American Statistical Association}, 80, 580-598.}
\keyword{datasets}
