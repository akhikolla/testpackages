\name{family}
\alias{family}

\title{Family of distributions}

\description{
Families of distributions supported by package \code{stepR}.

\bold{Deprecation warning:} This overviw is deprecated, but still given and up to date for some older, deprecated functions, however, may be removed in a future version. For an overview about the parametric families supported by the new functions see \code{\link{parametricFamily}}.
}

\details{
Package \code{stepR} supports several families of distributions (mainly exponential) to model the data, some of which require additional (fixed) parameters. In particular, the following families are available:

\describe{
\item{\code{"gauss"}}{normal distribution with unknown mean but known, fixed standard deviation given as a single \code{\link{numeric}} (will be estimated using \code{\link{sdrobnorm}} if omitted); cf. \code{\link{dnorm}}.}
\item{\code{"gaussvar"}}{normal distribution with unknown variance but known, fixed mean assumed to be zero; cf. \code{\link{dnorm}}.}
\item{\code{"poisson"}}{Poisson distribution with unknown intensity (no additional parameter); cf. \code{\link{dpois}}.}
\item{\code{"binomial"}}{binomial distribution with unknown success probability but known, fixed size given as a single \code{\link{integer}}; cf. \code{\link{dbinom}}.}
\item{\code{"gaussKern"}}{normal distribution with unknown mean and unknown, fixed standard deviation (being estimated using \code{\link{sdrobnorm}}), after filtering with a fixed filter which needs to be given as the additional parameter (a \code{\link{dfilter}} object); cf. \code{\link{dfilter}}.}
}

The family is selected via the \code{family} argument, providing the corresponding string, while the \code{param} argument contains the parameters if any.
}

\note{Beware that not all families can be chosen for all functions.}

\examples{
# illustrating different families fitted to the same binomial data set
size <- 200
n <- 200
# truth
p <- 10^seq(-3, -0.1, length = n)
# data
y <- rbinom(n, size, p)
plot(y)
lines(size * p, col = "red")
# fit 4 jumps, binomial family
jumps <- 4
bfit <- steppath(y, family = "binomial", param = size, max.blocks = jumps)
lines(bfit[[jumps]], col = "orange")
# Gaussian approximation with estimated variance
gfit <- steppath(y, family = "gauss", max.blocks = jumps)
lines(gfit[[jumps]], col = "green3", lty = 2)
# Poisson approximation
pfit <- steppath(y, family = "poisson", max.blocks = jumps)
lines(pfit[[jumps]], col = "blue", lty = 2)
legend("topleft", legend = c("binomial", "gauss", "poisson"), lwd = 2,
  col = c("orange", "green3", "blue"))
}

\seealso{\link{Distributions}, \code{\link{parametricFamily}}, \code{\link{dnorm}}, \code{\link{dpois}}, \code{\link{dbinom}}, \code{\link{dfilter}}, \code{\link{sdrobnorm}}}

\keyword{distribution}
