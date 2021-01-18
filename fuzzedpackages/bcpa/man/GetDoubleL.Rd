\name{GetDoubleL}
\alias{GetDoubleL}
\title{Obtain log-likelihood and parameter estimates for a given break point.}
\usage{
GetDoubleL(x, t, tbreak, ...)
}
\arguments{
  \item{x}{vector of time series values.}

  \item{t}{vector of times of measurements associated with
  x.}

  \item{tbreak}{breakpoint to test (in terms of the INDEX
  within "t" and "x", not actual time value).}

  \item{...}{additional parameters to pass to
  \code{\link{GetRho}}.}
}
\value{
a vector containing the parameters and the negative
log-likelihoods in order: \code{mu1, sigma1, tau1, LL1,
mu2, sigma2, tau2, LL2}
}
\description{
Takes a time series with values \code{x} obtained at time
\code{t} and a time break \code{tbreak}, and returns the
estimates of \eqn{\mu}, \eqn{\sigma} and \eqn{\tau} (or
\eqn{\rho}) as well as the negative log-likelihood of those
estimates before and after the break. Mostly for use
internally within \code{\link{GetBestBreak}}.
}
\author{
Eliezer Gurarie
}
\seealso{
\code{\link{GetBestBreak}} uses this function, while this
function uses \code{\link{GetRho}}
}

