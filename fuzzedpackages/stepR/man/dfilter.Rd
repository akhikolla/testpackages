\name{dfilter}
\alias{dfilter}
\alias{print.dfilter}

\title{Digital filters}

\description{
Create digital filters.

\bold{Deprecation warning:} This function is mainly used for patchlamp recordings and may be transferred to a specialised package.
}

\usage{
dfilter(type = c("bessel", "gauss", "custom"), param = list(pole = 4, cutoff = 1 / 10),
  len = ceiling(3/param$cutoff))
\method{print}{dfilter}(x, \dots)
}

\arguments{
\item{type}{allows to choose Bessel, Gauss or custom filters}
\item{param}{for a \code{"bessel"} filter a \code{\link{list}} with entries \code{pole} and \code{cutoff} giving the filter's number of poles (order) and cut-off frequency, resp.; for a \code{"gauss"} filter the filter's bandwidth (standard deviation) as a single \code{\link{numeric}}; for a custom filter either a \code{\link{numeric}} vector specifying the filter's kernel or a list with items \code{kern} and \code{step} of the same length giving the filter's kernel and step-response, resp.}
\item{len}{filter length (unnecessary for \code{"custom"} filters}
\item{x}{the object}
\item{\dots}{for generic methods only}
}

\value{
Returns a list of \code{\link{class}} \code{dfilter} that contains elements \code{kern} and \code{step}, the (digitised) filter kernel and step-response, resp., as well as an element \code{param} containing the argument \code{param}, for a \code{"bessel"} filter alongside the corresponding analogue kernel, step response, power spectrum, and autocorrelation function depending on time or frequency as elements \code{kernfun}, \code{stepfun}, \code{spectrum}, and \code{acfun}, resp.
}

\seealso{\code{\link{filter}}, \code{\link{convolve}}, \code{\link{BesselPolynomial}}, \code{\link{Normal}}, \link{family}}

\examples{
# 6-pole Bessel filter with cut-off frequency 1 / 100, with length 100 (too short!)
dfilter("bessel", list(pole = 6, cutoff = 1 / 100), 100)
# custom filter: running mean of length 3
dfilter("custom", rep(1, 3))
dfilter("custom", rep(1, 3))$kern # normalised!
dfilter("custom", rep(1, 3))$step
# Gaussian filter with bandwidth 3 and length 11 (from -5 to 5)
dfilter("gauss", 3, 11)
}

\keyword{ts}
