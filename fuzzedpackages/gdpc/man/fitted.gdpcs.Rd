\name{fitted.gdpcs}
\alias{fitted.gdpcs}
\title{Get Reconstructed Time Series From a gdpcs Object}
\description{
  Get reconstructed time series from a gdpcs object.
}


\usage{
  \method{fitted}{gdpcs}(object, num_comp = 1, \dots)
}
\arguments{
  \item{object}{An object of class \code{gdpcs}, usually the result of \link{auto.gdpc}.}
  \item{num_comp}{Integer indicating how many components to use for the reconstruction. Default is 1.}
  \item{\dots}{Additional arguments for compatibility.}
}


\value{ A matrix that is the reconstruction of the original series.}
\author{
  Daniel Peña, Ezequiel Smucler, Victor Yohai
}



\seealso{
  \code{\link{gdpc}}, \code{\link{auto.gdpc}}, \code{\link{plot.gdpc}}
}

\examples{
  T <- 200 #length of series
  m <- 200 #number of series
  set.seed(1234)
  f <- rnorm(T + 1)
  x <- matrix(0, T, m)
  u <- matrix(rnorm(T * m), T, m)
  for (i in 1:m) {
    x[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
  }
  #Choose number of lags using the LOO criterion.
  #k_max=2 to keep computation time low
  autofit <- auto.gdpc(x, k_max = 2, auto_comp = FALSE, num_comp = 2) 
  recons <- fitted(autofit, num_comp = 2)
}

\keyword{ts}