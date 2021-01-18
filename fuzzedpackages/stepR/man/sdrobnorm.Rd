\name{sdrobnorm}
\alias{sdrobnorm}

\title{Robust standard deviation estimate}
\description{Robust estimation of the standard deviation of Gaussian data.}

\usage{
sdrobnorm(x, p = c(0.25, 0.75), lag = 1,
          supressWarningNA = FALSE, supressWarningResultNA = FALSE)
}

\arguments{
\item{x}{a vector of numerical observations. \code{NA} entries will be removed with a warning. The warning can be supressed by setting \code{supressWarningNA} to \code{TRUE}. Other non finite values are not allowed}
\item{p}{vector of two distinct probabilities}
\item{lag}{a single integer giving the lag of the difference used, see \code{\link{diff}}, if a numeric is passed a small tolerance will be added and the value will be converted by \code{\link{as.integer}}}
\item{supressWarningNA}{a single logical, if \code{TRUE} no warning will be given for \code{NA} entries in \code{x}}
\item{supressWarningResultNA}{a single logical, if \code{TRUE} no warning will be given if the result is \code{NA}}
}

\value{
Returns the estimate of the sample's standard deviation, i.e. a single non-negative numeric, \code{NA} if \code{length(x) < lag + 2}.
}

\details{
Compares the difference between the estimated sample quantile corresponding to \code{p} after taking (\code{lag}ged) differences) with the corresponding theoretical quantiles of Gaussian white noise to determine the standard deviation under a Gaussian assumption. If the data contain (few) jumps, this will (on average) be a slight overestimate of the true standard deviation.

This estimator has been inspired by (1.7) in (\cite{Davies and Kovac}, 2001).
}



\seealso{\code{\link{sd}}, \code{\link{diff}}, \link{parametricFamily}, \link{family}}

\references{
Davies, P. L., Kovac, A. (2001) Local extremes, runs, strings and multiresolution. \emph{The Annals of Statistics} \bold{29}, 1--65.
}

\examples{
# simulate data sample
y <- rnorm(100, c(rep(1, 50), rep(10, 50)), 2)
# estimate standard deviation
sdrobnorm(y)
}

\keyword{nonparametric}
