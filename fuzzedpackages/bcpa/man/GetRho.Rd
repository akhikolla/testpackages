\name{GetRho}
\alias{GetRho}
\alias{GetRho2}
\title{Characteristic time / auto-correlation for irregular time series}
\usage{
GetRho(x, t, tau = TRUE)
}
\arguments{
  \item{x}{vector of time series values.}

  \item{t}{vector of times of measurements associated with
  x.}

  \item{tau}{whether or not to estimate time scale
  \eqn{\tau} (preferred) or autocorrelation \eqn{\rho}.}
}
\value{
Returns a vector of length two: the estimate and the
negative log-likelihood
}
\description{
Estimates characteristic time \eqn{\tau} or
auto-correlation \eqn{\rho} from a gappy time series
dataset.  It works first by estimating the mean and
standard deviation directly from the time series X, using
these to standardize the time series, and then optimizes
for the likelihood of the value of \eqn{\tau} or
\eqn{\rho}.
}
\author{
Eliezer Gurarie
}

