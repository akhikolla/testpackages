\name{stepcand}
\alias{stepcand}

\title{Forward selection of candidate jumps}

\description{
Find candidates for jumps in serial data by forward selection.
}

\usage{
stepcand(y, x = 1:length(y), x0 = 2 * x[1] - x[2], max.cand = NULL,
  family = c("gauss", "gaussvar", "poisson", "binomial", "gaussKern"), param = NULL,
  weights = rep(1, length(y)), cand.radius = 0)
}

\arguments{
\item{y}{a numeric vector containing the serial data}
\item{x}{a numeric vector of the same length as \code{y} containing the corresponding sample points}
\item{x0}{a single numeric giving the last unobserved sample point directly before sampling started}
\item{max.cand}{single integer giving the maximal number of blocks to find; defaults to using all data (note: there will be one block more than the number of jumps}
\item{family}{distribution of the errors, either \code{"gauss"}, \code{"poisson"} or \code{"binomial"}; \code{"gaussInhibit"} is like \code{"gauss"} forbids jumps getting close together or to the ends in \code{\link{steppath.stepcand}}, "gaussInhibitBoth" already forbids this in \code{stepcand} (not recommended)}
\item{param}{additional parameters specifying the distribution of the errors; the number of trials for family \code{"binomial"}; for \code{gaussInhibit} and \code{gaussInhibitBoth} a numeric of length 3 with components \code{"start"}, \code{"middle"} and \code{"end"} preventing the first jump from getting closer to \code{x0} than the \code{"start"} value, any two jumps from getting closer than the \code{"middle"} value, and the last jump from geting closer than the \code{"end"} value to the end, all distances measured by weights (cf. example below)}
\item{weights}{a numeric vector of the same length as \code{y} containing non-negative weights}
\item{cand.radius}{a non-negative integer: adds for each candidate found all indices that are at most \code{cand.radius} away}
}

\value{
An object of class \code{stepcand} extending class \code{\link{stepfit}} such that it can be used as an input to \code{\link{steppath.stepcand}}: additionally contains columns
\item{\code{cumSum}}{The cumulative sum of \code{x} up to \code{rightEnd}.}
\item{\code{cumSumSq}}{The cumulative sum of squares of \code{x} up to \code{rightEnd} (for \code{family = "gauss"}).}
\item{\code{cumSumWe}}{The cumulative sum of weights up to \code{rightEnd}.}
\item{\code{improve}}{The improvement this jump brought about when it was selected.}
}

\seealso{\code{\link{steppath}}, \code{\link{stepfit}}, \link{family}}

\examples{
# simulate 5 blocks (4 jumps) within a total of 100 data points
b <- c(sort(sample(1:99, 4)), 100)
f <- rep(rnorm(5, 0, 4), c(b[1], diff(b)))
rbind(b = b, f = unique(f), lambda = exp(unique(f) / 10) * 20)
# add gaussian noise
x <- f + rnorm(100)
# find 10 candidate jumps
stepcand(x, max.cand = 10)
# for poisson observations
y <- rpois(100, exp(f / 10) * 20)
# find 10 candidate jumps
stepcand(y, max.cand = 10, family = "poisson")
# for binomial observations
size <- 10
z <- rbinom(100, size, pnorm(f / 10))
# find 10 candidate jumps
stepcand(z, max.cand = 10, family = "binomial", param = size)
}

\keyword{nonparametric}
