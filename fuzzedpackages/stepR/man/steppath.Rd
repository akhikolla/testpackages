\name{steppath}
\alias{steppath}
\alias{steppath.default}
\alias{steppath.stepcand}
\alias{[[.steppath}
\alias{length.steppath}
\alias{print.steppath}
\alias{logLik.steppath}

\title{Solution path of step-functions}

\description{
Find optimal fits with step-functions having jumps at given candidate positions for all possible subset sizes.
}

\usage{
steppath(y, \dots, max.blocks)
\method{steppath}{default}(y, x = 1:length(y), x0 = 2 * x[1] - x[2], max.cand = NULL,
  family = c("gauss", "gaussvar", "poisson", "binomial", "gaussKern"), param = NULL,
  weights = rep(1, length(y)), cand.radius = 0, \dots, max.blocks = max.cand)
\method{steppath}{stepcand}(y, \dots, max.blocks = sum(!is.na(y$number)))
\method{[[}{steppath}(x, i)
\method{length}{steppath}(x)
\method{print}{steppath}(x, \dots)
\method{logLik}{steppath}(object, df = NULL, nobs = object$cand$rightIndex[nrow(object$cand)], \dots)
}

\arguments{
for \code{steppath}:
\item{y}{either an object of class \code{\link{stepcand}} for \code{steppath.stepcand} or a numeric vector containing the serial data for \code{steppath.default}}
\item{x, x0, max.cand, family, param, weights, cand.radius}{for \code{steppath.default} which calls \code{\link{stepcand}}; see there}
\item{max.blocks}{single integer giving the maximal number of blocks to find; defaults to number of candidates (note: there will be one block more than the number of jumps}
\item{\dots}{for generic methods only}
for methods on a \code{steppath} object \code{x} or \code{object}:
\item{object}{the object}
\item{i}{if this is an integer returns the fit with \code{i} blocks as an object of class \code{\link{stepcand}}, else the standard behaviour of a \code{\link{list}}}
\item{df}{the number of estimated parameters: by default the number of blocks for families \code{poisson} and \code{binomial}, one more (for the variance) for family \code{gauss}}
\item{nobs}{the number of observations used for estimating}
}

\value{
For \code{steppath} an object of class \code{steppath}, i.e. a \code{\link{list}} with components
\item{\code{path}}{A list of length \code{length(object)} where the \code{i}th element contains the best fit by a step-function having \code{i-1} jumps (i.e. \code{i} blocks), given by the candidates indices}
\item{\code{cost}}{A numeric vector of length \code{length(object)} giving the value of the cost functional corresponding to the solutions.}
\item{\code{cand}}{An object of class \code{\link{stepcand}} giving the candidates among which the jumps were selected.}
\code{[[.steppath} returns the fit with \code{i} blocks as an object of class \code{\link{stepfit}}; \code{length.steppath} the maximum number of blocks for which a fit has been computed. \code{logLik.stepfit} returns an object of class \code{\link{logLik}} giving the likelihood of the data given the fits corresponding to \code{cost}, e.g. for use with \code{\link{AIC}.}
}

\seealso{\code{\link{stepcand}}, \code{\link{stepfit}}, \link{family}, \code{\link{logLik}}, \code{\link{AIC}}}

\references{
Friedrich, F., Kempe, A., Liebscher, V., Winkler, G. (2008) Complexity penalized M-estimation: fast computation. \emph{Journal of Computational and Graphical Statistics} \bold{17}(1), 201--224.
}

\examples{
# simulate 5 blocks (4 jumps) within a total of 100 data points
b <- c(sort(sample(1:99, 4)), 100)
f <- rep(rnorm(5, 0, 4), c(b[1], diff(b)))
# add Gaussian noise
x <- f + rnorm(100)
# find 10 candidate jumps
cand <- stepcand(x, max.cand = 10)
cand
# compute solution path
path <- steppath(cand)
path
plot(x)
lines(path[[5]], col = "red")
# compare result having 5 blocks with truth
fit <- path[[5]]
fit
logLik(fit)
AIC(logLik(fit))
cbind(fit, trueRightEnd = b, trueLevel = unique(f))
# for poisson observations
y <- rpois(100, exp(f / 10) * 20)
# compute solution path, compare result having 5 blocks with truth
cbind(steppath(y, max.cand = 10, family = "poisson")[[5]],
  trueRightEnd = b, trueIntensity = exp(unique(f) / 10) * 20)
# for binomial observations
size <- 10
z <- rbinom(100, size, pnorm(f / 10))
# compute solution path, compare result having 5 blocks with truth
cbind(steppath(z, max.cand = 10, family = "binomial", param = size)[[5]],
  trueRightEnd = b, trueIntensity = pnorm(unique(f) / 10))
# an example where stepcand is not optimal but indices found are close to optimal ones
blocks <- c(rep(0, 9), 1, 3, rep(1, 9))
blocks
stepcand(blocks, max.cand = 3)[,c("rightEnd", "value", "number")]
# erroneously puts the "1" into the right block in the first step
steppath(blocks)[[3]][,c("rightEnd", "value")]
# putting the "1" in the middle block is optimal
steppath(blocks, max.cand = 3, cand.radius = 1)[[3]][,c("rightEnd", "value")]
# also looking in the 1-neighbourhood remedies the problem
}

\keyword{nonparametric}
